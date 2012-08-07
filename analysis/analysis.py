"""
Set of routines to analyse single-subject fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import os, os.path

import numpy as np

import fmri_tools.utils

import glass_coherence_block.config
import glass_coherence_block.analysis.paths
import glass_coherence_block.analysis.analysis

def exp_glm( paths, conf ):
	"""Experiment GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	stim_files = [ "%s%d.txt" % ( paths[ "ana" ][ "exp_time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = [ "%.2f" % coh for coh in conf[ "stim" ][ "coh_levels" ] ]

	trend_coef = [ [ -3, -1, +1, +3 ],  # linear
	               [ +1, -1, -1, +1 ],  # quadratic
	               [ -1, +3, -3, +1 ],  # cubic
	               [ +1, +1, +1, +1 ]   # mean
	             ]

	trend_lbl = [ "lin", "quad", "cub", "mean" ]

	trend_con = [ "SYM: " +
	              " ".join( [ "%d*%s" % con
	                          for con in zip( t_coef, stim_labels )
	                        ]
	                      )
	              for t_coef in trend_coef
	            ]

	for hemi in [ "lh", "rh" ]:

		fit_file = "%s_%s" % ( paths[ "ana" ][ "exp_fits" ], hemi )
		glm_file = "%s_%s" % ( paths[ "ana" ][ "exp_glm" ], hemi )
		beta_file = "%s_%s" % ( paths[ "ana" ][ "exp_beta" ], hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in paths[ "func" ][ "surf_files" ]
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", conf[ "ana" ][ "poly_ord" ],
		                  "-local_times",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-fitts", "%s.niml.dset" % fit_file,
		                  "-bucket", "%s.niml.dset" % glm_file,
		                  "-cbucket", "%s.niml.dset" % beta_file,
		                  "-jobs", "16",
		                  "-tout",
		                  "-overwrite",
		                  "-x1D_stop",
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

		for i_con in xrange( len( trend_con ) ):

			glm_cmd.extend( [ "-gltsym",
			                  trend_con[ i_con ]
			                ]
			              )

			glm_cmd.extend( [ "-glt_label",
			                  "%d" % ( i_con + 1 ),
			                  trend_lbl[ i_con ]
			                ]
			              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "%s.REML_cmd" % glm_file )

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

		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

	os.chdir( start_dir )


def loc_mask( paths, conf ):
	"""a"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	loc_stat_brik = "16"

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "exp_glm" ], hemi )

		loc_fdr_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_loc_fdr" ], hemi )

		# convert the statistics for the localiser to a q (FDR) value
		fdr_cmd = [ "3dFDR",
		            "-input", "%s[%s]" % ( glm_file, loc_stat_brik ),
		            "-prefix", loc_fdr_file,
		            "-qval",
		            "-float",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		loc_mask_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_loc_mask" ], hemi )

		q = conf[ "ana" ][ "loc_q" ]

		mask_cmd = [ "3dcalc",
		             "-a", loc_fdr_file,
		             "-b", "%s[%s]" % ( glm_file, loc_stat_brik ),
		             "-expr", "within( a, 0, %.6f ) * ispositive( b )" % q,
		             "-prefix", loc_mask_file,
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( mask_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_mask_file = "%s_%s-full" % ( paths[ "ana" ][ "exp_loc_mask" ], hemi )

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
	beta_briks = "48,49,50,51"

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	for hemi in [ "lh", "rh" ]:

		# dataset holding the beta weights
		beta_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "exp_beta" ], hemi )

		# design matrix file
		mat_file = os.path.join( paths[ "ana" ][ "exp_dir" ],
		                         "exp_design.xmat.1D"
		                       )

		# baseline timecourse dataset, to write
		bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_bltc" ], hemi )

		# generate an average baseline timecourse
		bl_cmd = [ "3dSynthesize",
		           "-cbucket", beta_file,
		           "-matrix", mat_file,
		           "-select", "baseline",
		           "-prefix", bltc_file,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( bl_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# baseline (point-estimate) dataset, to write
		bl_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_bl" ], hemi )

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
		psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_psc" ], hemi )

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
		full_psc_file = "%s_%s-full" % ( paths[ "ana" ][ "exp_psc" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( psc_file,
		                                 full_psc_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )


def roi_xtr( paths, conf ):
	"""Extract PSC and statistics data from ROIs"""

	for hemi in [ "lh", "rh" ]:

		loc_mask_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "exp_loc_mask" ],
		                                           hemi
		                                         )

		cmask = "-a %s -expr step(a)" % loc_mask_file

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "roi_dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			roi_psc_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			if os.path.exists( roi_psc_file ):
				os.remove( roi_psc_file )

			# our input dataset
			data_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "exp_psc" ],
			                                       hemi
			                                     )

			# use the ROI file to mask the input dataset
			xtr_cmd = [ "3dmaskdump",
			            "-mask", roi_file,
			            "-cmask", cmask,
			            "-mrange", roi_val, roi_val,
			            "-noijk",
			            "-o", roi_psc_file,
			            data_file
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


def group_rois( paths, conf ):
	"""a"""

	roi_names = [ roi_info[ 0 ] for roi_info in conf[ "ana" ][ "rois" ] ]

	for roi_name in roi_names:

		roi_path = "%s-%s.txt" % ( paths[ "roi_mean" ], roi_name )

		roi_file = open( roi_path, "w+" )

		for subj_id in conf[ "all_subj" ]:

			subj_conf = glass_coherence_block.config.get_conf( subj_id )

			subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

			roi_data = []

			for hemi in [ "lh", "rh" ]:

				roi_psc_file = "%s_%s_%s.txt" % ( subj_paths[ "rois" ][ "psc" ],
				                                  roi_name,
				                                  hemi
				                                )

				roi_hemi_data = np.loadtxt( roi_psc_file )

				roi_data.append( roi_hemi_data )

			roi_data = np.vstack( roi_data )

			roi_data = np.mean( roi_data, axis = 0 )

			subj_roi_txt = ( subj_id +
			                 "\t%.18e\t%.18e\t%.18e\t%.18e\n" % tuple( roi_data )
			               )

			roi_file.write( subj_roi_txt )

		roi_file.close()

