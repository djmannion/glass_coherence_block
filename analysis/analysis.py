"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
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
	             [ -1, +3, -3, +1 ],  # cubic
	             [ +1, +1, +1, +1 ]   # all
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
		                  for surf_file in paths[ "func" ][ "surf_files" ]
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
		                  "-gltsym", con[ 3 ],
		                  "-glt_label", "4", "all",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-jobs", "16",
		                  "-fitts", "%s.niml.dset" % fit_file,
		                  "-bucket", "%s.niml.dset" % glm_file,
		                  "-cbucket", "%s.niml.dset" % beta_file,
		                  "-tout",
		                  "-overwrite",
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		# convert the output to full
		for out_file in [ fit_file, glm_file, beta_file ]:

			fmri_tools.utils.sparse_to_full( "%s.niml.dset" % out_file,
			                                 "%s-full" % out_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )

	os.chdir( start_dir )
