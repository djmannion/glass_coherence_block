"""
Set of routines to analyse single-subject fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import os

import numpy as np

import fmri_tools.utils


def exp_glm( paths, conf ):
	"""Experiment GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	stim_files = [ "%s%d.txt" % ( paths[ "ana" ][ "exp_time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = [ "%.2f" % coh for coh in conf[ "stim" ][ "coh_levels" ] ]

	con_coef = [ [ -3, -1, +1, +3 ],  # linear
	             [ +1, -1, -1, +1 ],  # quadratic
	             [ -1, +3, -3, +1 ]  # cubic
	           ]


	con = [ "SYM: %d*%s %d*%s %d*%s %d*%s" % (
	          coef[ 0 ], stim_labels[ 0 ],
	          coef[ 1 ], stim_labels[ 1 ],
	          coef[ 2 ], stim_labels[ 2 ],
	          coef[ 3 ], stim_labels[ 3 ]
	        )
	        for coef in con_coef
	      ]

	for hemi in [ "lh", "rh" ]:

		fit_file = "%s_%s" % ( paths[ "ana" ][ "exp_fits" ], hemi )
		glm_file = "%s_%s" % ( paths[ "ana" ][ "exp_glm" ], hemi )
		beta_file = "%s_%s" % ( paths[ "ana" ][ "exp_beta" ], hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in paths[ "func" ][ "smooth_files" ]
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "4",
		                  "-num_stimts", "4",
		                  "-stim_label", "1", stim_labels[ 0 ],
		                  "-stim_times", "1", stim_files[ 0 ], "SPMG1(16)",
		                  "-stim_label", "2", stim_labels[ 1 ],
		                  "-stim_times", "2", stim_files[ 1 ], "SPMG1(16)",
		                  "-stim_label", "3", stim_labels[ 2 ],
		                  "-stim_times", "3", stim_files[ 2 ], "SPMG1(16)",
		                  "-stim_label", "4", stim_labels[ 3 ],
		                  "-stim_times", "4", stim_files[ 3 ], "SPMG1(16)",
		                  "-local_times",
		                  "-gltsym", con[ 0 ],
		                  "-glt_label", "1", "linear",
		                  "-gltsym", con[ 1 ],
		                  "-glt_label", "2", "quadratic",
		                  "-gltsym", con[ 2 ],
		                  "-glt_label", "3", "cubic",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-jobs", "16",
		                  "-fitts", "%s.niml.dset" % fit_file,
		                  "-bucket", "%s.niml.dset" % glm_file,
		                  "-cbucket", "%s.niml.dset" % beta_file,
		                  "-tout",
		                  "-overwrite",
		                  "-x1D_stop"
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-input", " ".join( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                                   for surf_file in paths[ "func" ][ "smooth_files" ]
		                                 ]
		                               ),
		             "-Rbeta", "%s_reml.niml.dset" % beta_file,
		             "-tout",
		             "-Rbuck", "%s_reml.niml.dset" % glm_file,
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_glm_file = "%s_reml-full" % ( glm_file, hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		pad_node = "ld141"

		fmri_tools.utils.sparse_to_full( "%s_reml.niml.dset" % glm_file,
		                                 full_glm_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def mask_nodes( paths, conf ):
	"""a"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "exp_glm" ], hemi )

		fdr_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_fdr" ], hemi )

		# first, need to convert the statistics to a q (FDR) value
		fdr_cmd = [ "3dFDR",
		            "-input", glm_file,
		            "-prefix", fdr_file,
		            "-qval",
		            "-float",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		p = conf[ "ana" ][ "loc_p" ]

		q_briks = [ 2, 4, 6, 8 ]

		fdr_data_files = [ "%s[%d]" % ( fdr_file, q_brik )
		                   for q_brik in q_briks
		                 ]

		t_briks = [ 2, 4, 6, 8 ]

		glm_data_files = [ "%s[%d]" % ( glm_file, t_brik )
		                   for t_brik in t_briks
		                 ]

		loc_mask_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_loc_mask" ],
		                                      hemi
		                                    )

		expr = ( "or( " +
		         "and( within( a, 0, %.6f ), posval( e ) )," % p +
		         "and( within( b, 0, %.6f ), posval( f ) )," % p +
		         "and( within( c, 0, %.6f ), posval( g ) )," % p +
		         "and( within( d, 0, %.6f ), posval( h ) )" % p +
		         ")"
		       )


		# then, need to apply conjunction test
		conj_cmd = [ "3dcalc",
		             "-a", fdr_data_files[ 0 ],
		             "-b", fdr_data_files[ 1 ],
		             "-c", fdr_data_files[ 2 ],
		             "-d", fdr_data_files[ 3 ],
		             "-e", glm_data_files[ 0 ],
		             "-f", glm_data_files[ 1 ],
		             "-g", glm_data_files[ 2 ],
		             "-h", glm_data_files[ 3 ],
		             "-expr", expr,
		             "-prefix", loc_mask_file,
		             "-overwrite"
		            ]

		fmri_tools.utils.run_cmd( conj_cmd,
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
	beta_briks = "60,61,62,63"

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
		            "-expr", "100 * b/a",
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
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "roi_dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			roi_psc_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "exp_roi_psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			# our input dataset - either psc or glm
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

