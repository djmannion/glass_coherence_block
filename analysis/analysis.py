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

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	stim_files = [ "%s%d.txt" % ( paths[ "ana" ][ "exp_time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = [ "%.2f" % coh for coh in conf[ "stim" ][ "coh_levels" ] ]

	trend_coef = [ [ -3, -1, +1, +3 ],  # linear
	               [ +1, -1, -1, +1 ],  # quadratic
	               [ -1, +3, -3, +1 ]   # cubic
	             ]

	trend_lbl = [ "lin", "quad", "cub" ]

	trend_con = [ "SYM: " +
	              " ".join( [ "%d*%s" % con
	                          for con in zip( t_coef, stim_labels )
	                        ]
	                      )
	              for t_coef in trend_coef
	            ]

	f_con = "SYM: +1*%s | +1*%s | +1*%s | +1*%s " % tuple( stim_labels )

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
		                  "-fout",
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

		# f contrast
		glm_cmd.extend( [ "-gltsym",
		                  f_con
		                ]
		              )

		glm_cmd.extend( [ "-glt_label",
		                  "%d" % ( len( trend_con ) + 1 ),
		                  "fAll"
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
		             "-fout",
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

	loc_F_brik = "30"

	loc_t_briks = [ "23", "25", "27", "29" ]

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "exp_glm" ], hemi )

		loc_fdr_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_loc_fdr" ], hemi )

		# convert the statistics for the localiser to a q (FDR) value
		fdr_cmd = [ "3dFDR",
		            "-input", "%s[%s]" % ( glm_file, loc_F_brik ),
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

		# has to have a significant F statistic AND have at least one of it's
		# constituent t-values be positive
		expr = """within( a, 0, %.6f ) * or( ispositive( b ), ispositive( c ),
ispositive( d ), ispositive( e ) )""" % q

		mask_cmd = [ "3dcalc",
		             "-a", loc_fdr_file,
		             "-b", "%s[%s]" % ( glm_file, loc_t_briks[ 0 ] ),
		             "-c", "%s[%s]" % ( glm_file, loc_t_briks[ 1 ] ),
		             "-d", "%s[%s]" % ( glm_file, loc_t_briks[ 2 ] ),
		             "-e", "%s[%s]" % ( glm_file, loc_t_briks[ 3 ] ),
		             "-expr", expr,
		             "-prefix", loc_mask_file,
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( mask_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

	os.chdir( start_dir )


def trends( paths, conf ):
	"""a"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	loc_stat_brik = "16"

	trend_stat_briks = "10,12,14"

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "exp_glm" ], hemi )

		loc_mask_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_loc_mask" ], hemi )

		trend_fdr_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_trend_fdr" ], hemi )

		fdr_cmd = [ "3dFDR",
		            "-input", "%s[%s]" % ( glm_file, trend_stat_briks ),
		            "-prefix", trend_fdr_file,
		            "-mask", loc_mask_file,
		            "-float",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		 



def beta_to_psc( paths, conf ):
	"""Convert the GLM beta weights into units of percent signal change"""

	# these are the indices into the beta files for the data we want to convert
	beta_briks = "60,61,62,63"

	beta_briks = "9,15"

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

		avg_cmd = [ "3dMean",
		            "-prefix", bl_file,
		            "%s[0,5,10,15,20,25,30,35,40,45,50,55]" % beta_file
		          ]

		fmri_tools.utils.run_cmd( avg_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# dataset to hold the percent signal change, to write
		psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_psc" ], hemi )
		beta_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "exp_glm" ], hemi )

		# the input beta file, with sub-brick selector
		beta_sel = "%s[%s]" % ( beta_file, beta_briks )

		# check that the label is as expected
#		beta_label = fmri_tools.utils.get_dset_label( beta_sel )

#		assert( beta_label == [ "0.00#0", "0.33#0", "0.66#0", "1.00#0" ] )

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

		pad_node = "ld141"

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

