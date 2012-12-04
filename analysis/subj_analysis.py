"""
Set of routines to analyse single-subject fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import os, os.path
import logging

import numpy as np

import fmri_tools.utils


def glm( conf, paths ):
	"""Experiment GLM"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running GLM..." )

	start_dir = os.getcwd()

	os.chdir( paths.ana.base.full() )

	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	# label for each condition
	stim_labels = [ "{coh:.0f}".format( coh = ( coh * 100 ) )
	                for coh in conf[ "stim" ][ "coh_levels" ]
	              ]

	# files containing the stimulus times for each condition
	stim_times = [ paths.ana.stim_times.full( "_{c:d}.txt".format( c = cond_num ) )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	# contrast coefficients; each column corresponds to a condition
	con_coef = [ [ +1, +1, +1, +1 ],  # localiser
	             [ -3, -1, +1, +3 ]   # linear
	           ]

	# label for each contrast
	con_lbl = [ "stim", "lin" ]

	# put the contrasts into an AFNI-aware format
	con_str = [ "SYM: " +
	            " ".join( [ "{c:d}*{l:s}".format( c = c, l = l )
	                        for ( c, l ) in zip( coef, stim_labels )
	                      ]
	                    )
	            for coef in con_coef
	          ]

	# minus one because the range is inclusive
	censor_vols = conf[ "exp" ][ "pre_len_s" ] / conf[ "acq" ][ "tr_s" ] - 1

	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-{v:.0f}".format( v = censor_vols )

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ surf_file.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		                  for surf_file in paths.func.surfs
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf[ "acq" ][ "tr_s" ] ),
		                  "-polort", "{p:d}".format( p = conf[ "ana" ][ "poly_ord" ] ),
		                  "-ortvec", paths.summ.motion.full( ".txt" ), "mot",
		                  "-local_times",
		                  "-CENSORTR", censor_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "{n:d}".format( n = n_cond )
		                ]
		              )

		for i_stim in xrange( n_cond ):

			glm_cmd.extend( [ "-stim_label",
			                  "{sl:d}".format( sl = ( i_stim + 1 ) ),
			                  stim_labels[ i_stim ]
			                ]
			              )

			glm_cmd.extend( [ "-stim_times",
			                  "{sl:d}".format( sl = ( i_stim + 1 ) ),
			                  stim_times[ i_stim ],
			                  hrf_model
			                ]
			              )

		# loop through all our contrasts
		for i_con in xrange( len( con_coef ) ):

			glm_cmd.extend( [ "-gltsym",
			                  "'" + con_str[ i_con ] + "'"
			                ]
			              )

			glm_cmd.extend( [ "-glt_label",
			                  "{cl:d}".format( cl = ( i_con + 1 ) ),
			                  con_lbl[ i_con ]
			                ]
			              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", paths.ana.beta.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             "-tout",
		             "-Rbuck", paths.ana.glm.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( "'" +
		                 " ".join( [ surf_file.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		                             for surf_file in paths.func.surfs
		                           ]
		                         ) +
		                 "'"
		               )

		# run the proper GLM
		fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )


	os.chdir( start_dir )


def loc_mask( paths, conf ):
	"""Creates a mask from the localiser contrast"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	# brik in the glm file that contains the localiser statistic
	# this is verified below
	loc_stat_brik = "10"

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "glm" ], hemi )
		glm_file += "[%s]" % loc_stat_brik

		# check the localiser brik is as expected
		assert( fmri_tools.utils.get_dset_label( glm_file ) == [ "stim#0_Tstat" ] )

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


def loc_mask( conf, paths ):
	"""Form a mask from the GLM output"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running localising mask creation..." )

	# these are the indices into the beta files for the data we want to convert
	loc_brick = "[10]"

	for hemi in [ "lh", "rh" ]:

		glm_path = ( paths.ana.glm.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ) +
		             loc_brick
		           )

		# check that the label is as expected
		glm_label = fmri_tools.utils.get_dset_label( glm_path )
		assert( glm_label == [ "stim#0_Tstat" ] )

		fdr_path = paths.ana.fdr.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		mask_path = paths.ana.mask.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		fmri_tools.utils.loc_mask( glm_path = glm_path,
		                           fdr_path = fdr_path,
		                           mask_path = mask_path,
		                           q_thresh = conf[ "ana" ][ "loc_q" ],
		                           pos_only = True
		                         )

		pad_mask_path = paths.ana.mask.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )
		pad_nodes = "{nk:d}".format( nk = conf[ "subj" ][ "node_k" ][ hemi ] )

		fmri_tools.utils.sparse_to_full( in_dset = mask_path,
		                                 out_dset = pad_mask_path,
		                                 pad_node = pad_nodes
		                               )


def beta_to_psc( conf, psc ):
	"""Convert the GLM beta weights into units of percent signal change"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running beta to PSC conversion..." )

	# these are the indices into the beta files for the data we want to convert
	beta_bricks = "[60,61,62,63]"

	for hemi in [ "lh", "rh" ]:

		beta_path = paths.ana.beta.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		# check that the label is as expected
		beta_label = fmri_tools.utils.get_dset_label( beta_path + beta_bricks )
		assert( beta_label == [ "0#0", "33#0", "66#0", "100#0" ] )

		design_path = ( paths.ana.base + "exp_design.xmat.1D" ).full()

		bltc_path = paths.ana.bltc.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		bl_path = paths.ana.bl.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		psc_path = paths.ana.psc.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		fmri_tools.utils.beta_to_psc( beta_path = beta_path,
		                              beta_bricks = beta_bricks,
		                              design_path = design_path,
		                              bltc_path = bltc_path,
		                              bl_path = bl_path,
		                              psc_path = psc_path
		                            )

		pad_psc_path = paths.ana.psc.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )
		pad_nodes = "{nk:d}".format( nk = conf[ "subj" ][ "node_k" ][ hemi ] )

		fmri_tools.utils.sparse_to_full( in_dset = psc_path,
		                                 out_dset = pad_psc_path,
		                                 pad_node = pad_nodes
		                               )


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
