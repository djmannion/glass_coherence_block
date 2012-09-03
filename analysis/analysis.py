"""
Set of routines to analyse single-subject fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import os, os.path
import tempfile

import numpy as np
import scipy.stats

import fmri_tools.utils

import glass_coherence_block.config
import glass_coherence_block.analysis.paths


def glm( paths, conf ):
	"""Experiment GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	stim_files = [ "%s%d.txt" % ( paths[ "ana" ][ "time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = [ "%.2f" % coh for coh in conf[ "stim" ][ "coh_levels" ] ]

	# contrast coefficients; each column corresponds to a condition
	con_coef = [ [ +1, +1, +1, +1 ]
	           ]

	con_lbl = [ "stim" ]

	# put the contrast into an AFNI-aware format
	con_str = [ "SYM: " +
	            " ".join( [ "%d*%s" % con
	                        for con in zip( coef, stim_labels )
	                      ]
	                    )
	            for coef in con_coef
	          ]

	# minus one because the range is inclusive
	censor_vols = conf[ "exp" ][ "pre_len_s" ] / conf[ "acq" ][ "tr_s" ] - 1

	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-%d" % censor_vols

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in paths[ "func" ][ "surf_files" ]
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "-1",  # we pass our own below
		                  "-ortvec", paths[ "ana" ][ "bl_poly" ], "poly",
		                  "-ortvec", paths[ "summ" ][ "mot_est_file" ], "mot",
		                  "-local_times",
		                  "-CENSORTR", censor_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "%d" % n_cond
		                ]
		              )

		for i_stim in xrange( n_cond ):

			glm_cmd.extend( [ "-stim_label",
			                  "%d" % ( i_stim + 1 ),
			                  stim_labels[ i_stim ]
			                ]
			              )

			glm_cmd.extend( [ "-stim_times",
			                  "%d" % ( i_stim + 1 ),
			                  stim_files[ i_stim ],
			                  hrf_model
			                ]
			              )

		# loop through all our contrasts
		for i_con in xrange( len( con_coef ) ):

			glm_cmd.extend( [ "-gltsym",
			                  con_str[ i_con ]
			                ]
			              )

			glm_cmd.extend( [ "-glt_label",
			                  "%d" % ( i_con + 1 ),
			                  con_lbl[ i_con ]
			                ]
			              )

		# run this first GLM
		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		glm_file = "%s_%s" % ( paths[ "ana" ][ "glm" ], hemi )
		beta_file = "%s_%s" % ( paths[ "ana" ][ "beta" ], hemi )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", "%s_reml.niml.dset" % beta_file,
		             "-tout",
		             "-Rbuck", "%s_reml.niml.dset" % glm_file,
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( " ".join( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                             for surf_file in paths[ "func" ][ "surf_files" ]
		                           ]
		                         )
		               )

		# run the proper GLM
		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

	os.chdir( start_dir )


def loc_mask( paths, conf ):
	"""Creates a mask from the localiser contrast"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	# brik in the glm file that contains the localiser statistic
	# this is verified below
	loc_stat_brik = "16"

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "glm" ], hemi )
		glm_file += "[%s]" % loc_stat_brik

		# check the localiser brik is as expected
		assert( fmri_tools.utils.get_dset_label( glm_file ) == [ "mean#0_Tstat" ] )

		# to write
		loc_fdr_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_fdr" ], hemi )

		# convert the statistics for the localiser to a q (FDR) value
		fdr_cmd = [ "3dFDR",
		            "-input", glm_file,
		            "-prefix", loc_fdr_file,
		            "-qval",  # specify that we want q, not z
		            "-float",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		loc_mask_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_mask" ], hemi )

		q_thresh = conf[ "ana" ][ "loc_q" ]

		# create a localiser mask as nodes that both have a q that is below
		# threshold and have positive beta weights
		mask_cmd = [ "3dcalc",
		             "-a", loc_fdr_file,
		             "-b", glm_file,
		             "-expr", "within( a, 0, %.6f ) * ispositive( b )" % q_thresh,
		             "-prefix", loc_mask_file,
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( mask_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_mask_file = "%s_%s-full" % ( paths[ "ana" ][ "loc_mask" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( loc_mask_file,
		                                 full_mask_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def beta_to_psc( paths, conf ):
	"""Convert the GLM beta weights into units of percent signal change"""

	# these are the indices into the beta files for the data we want to convert
	beta_briks = "0,1,2,3"

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	for hemi in [ "lh", "rh" ]:

		# dataset holding the beta weights
		beta_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

		# design matrix file
		mat_file = os.path.join( paths[ "ana" ][ "base_dir" ],
		                         "exp_design.xmat.1D"
		                       )

		# baseline timecourse dataset, to write
		bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bltc" ], hemi )

		# generate an average baseline timecourse
		bl_cmd = [ "3dSynthesize",
		           "-cbucket", beta_file,
		           "-matrix", mat_file,
		           "-cenfill", "none",  # this is important
		           "-select", "baseline",
		           "-prefix", bltc_file,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( bl_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# baseline (point-estimate) dataset, to write
		bl_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bl" ], hemi )

		# average baseline timecourse across time
		avg_cmd = [ "3dTstat",
		            "-mean",
		            "-overwrite",
		            "-prefix", bl_file,
		            bltc_file
		          ]

		fmri_tools.utils.run_cmd( avg_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# dataset to hold the percent signal change, to write
		psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "psc" ], hemi )

		# the input beta file, with sub-brick selector
		beta_sel = "%s[%s]" % ( beta_file, beta_briks )

		# check that the label is as expected
		beta_label = fmri_tools.utils.get_dset_label( beta_sel )
		assert( beta_label == [ "0.00#0", "0.33#0", "0.66#0", "1.00#0" ] )

		# compute psc
		# from http://afni.nimh.nih.gov/sscc/gangc/TempNorm.html
		psc_cmd = [ "3dcalc",
		            "-fscale",
		            "-a", bl_file,
		            "-b", beta_sel,
		            "-expr", "100 * b/a * step (1- abs(b/a))",
		            "-prefix", psc_file,
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( psc_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_psc_file = "%s_%s-full" % ( paths[ "ana" ][ "psc" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( psc_file,
		                                 full_psc_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def roi_xtr( paths, conf ):
	"""Extract PSC and statistics data from ROIs"""

	for hemi in [ "lh", "rh" ]:

		# the *full* localiser mask file
		loc_mask_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "loc_mask" ],
		                                           hemi
		                                         )

		# expression to apply the localiser mask
		cmask_expr = "-a %s -expr step(a)" % loc_mask_file

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			roi_psc_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			# 3dmaskdump won't overwrite, so need to manually remove any prior file
			if os.path.exists( roi_psc_file ):
				os.remove( roi_psc_file )

			# our input dataset
			data_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "psc" ],
			                                       hemi
			                                     )

			# use the ROI file to mask the input dataset
			xtr_cmd = [ "3dmaskdump",
			            "-mask", roi_file,
			            "-cmask", cmask_expr,
			            "-mrange", roi_val, roi_val,
			            "-noijk",
			            "-o", roi_psc_file,
			            data_file
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


def raw_adj( paths, conf ):
	"""Concatenates raw timecourses and adjusts them for baselines"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	# since we use censoring in the GLM rather than trimming the data, we need to
	# specify a volume range here
	censor_vols = conf[ "exp" ][ "pre_len_s" ] / conf[ "acq" ][ "tr_s" ]
	brick_range = "[%d..$]" % censor_vols

	for hemi in [ "lh", "rh" ]:

		# create a raw input dataset, concatentated across all runs
		surf_files = [ "%s_%s.niml.dset%s" % ( surf_file, hemi, brick_range )
		               for surf_file in paths[ "func" ][ "surf_files" ]
		             ]

		raw_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "raw" ], hemi )

		cat_cmd = [ "3dTcat",
		            "-overwrite",
		            "-prefix", raw_file
		          ]

		cat_cmd.extend( surf_files )

		fmri_tools.utils.run_cmd( cat_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# adjust the raw file by subtracting the baseline
		raw_adj_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "raw_adj" ], hemi )

		bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bltc" ], hemi )

		adj_cmd = [ "3dcalc",
		            "-fscale",
		            "-a", bltc_file,
		            "-b", raw_file,
		            "-expr", "b-a",
		            "-prefix", raw_adj_file,
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( adj_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# dataset holding the beta weights
		beta_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

		# design matrix file
		mat_file = os.path.join( paths[ "ana" ][ "base_dir" ],
		                         "exp_design.xmat.1D"
		                       )

		pred_adj_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "pred_adj" ], hemi )

		# generate an signal timecourse
		bl_cmd = [ "3dSynthesize",
		           "-cbucket", beta_file,
		           "-matrix", mat_file,
		           "-cenfill", "none",  # important to match the censored data
		           "-select", "allfunc",
		           "-prefix", pred_adj_file,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( bl_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_raw_file = "%s_%s-full" % ( paths[ "ana" ][ "raw" ], hemi )
		full_raw_adj_file = "%s_%s-full" % ( paths[ "ana" ][ "raw_adj" ], hemi )
		full_pred_adj_file = "%s_%s-full" % ( paths[ "ana" ][ "pred_adj" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		pad_files = [ [ raw_file, full_raw_file ],
		              [ raw_adj_file, full_raw_adj_file ],
		              [ pred_adj_file, full_pred_adj_file ]
		            ]

		for ( sparse_file, full_file ) in pad_files:

			fmri_tools.utils.sparse_to_full( sparse_file,
			                                 full_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )

	os.chdir( start_dir )


def roi_tc( paths, conf ):
	"""Compile raw and predicted (adjusted) timecourses for each ROI."""

	for hemi in [ "lh", "rh" ]:

		# the *full* localiser mask file
		loc_mask_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "loc_mask" ],
		                                           hemi
		                                         )

		# expression to apply the localiser mask
		cmask_expr = "-a %s -expr step(a)" % loc_mask_file

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			roi_raw_adj_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "raw_adj_tc" ],
			                                      roi_name,
			                                      hemi
			                                    )

			raw_adj_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "raw_adj" ],
			                                          hemi
			                                        )

			roi_pred_adj_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "pred_adj_tc" ],
			                                       roi_name,
			                                       hemi
			                                     )

			pred_adj_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "pred_adj" ],
			                                           hemi
			                                         )

			data_files = [ [ roi_raw_adj_file, raw_adj_file ],
			               [ roi_pred_adj_file, pred_adj_file ]
			             ]

			for ( out_file, in_file ) in data_files:

				# 3dmaskdump won't overwrite, so need to manually remove any prior data
				if os.path.exists( out_file ):
					os.remove( out_file )

				# use the ROI file to mask the input dataset
				xtr_cmd = [ "3dmaskdump",
				            "-mask", roi_file,
				            "-cmask", cmask_expr,
				            "-mrange", roi_val, roi_val,
				            "-noijk",
				            "-o", out_file,
				            in_file
				          ]

				fmri_tools.utils.run_cmd( xtr_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )
