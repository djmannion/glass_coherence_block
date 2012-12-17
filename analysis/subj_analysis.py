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


def beta_to_psc( conf, paths ):
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


def roi_prep( conf, paths ):
	"""Prepares the ROI datasets"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running ROI preparation..." )

	start_dir = os.getcwd()

	os.chdir( paths.roi.base.dir() )

	for hemi in [ "lh", "rh" ]:

		# first step is to extract a subset of the visual localiser rois
		vl_path = paths.roi.vl.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )
		vl_sub_path = paths.roi.vl_subset.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		vl_expr = ( "'a*amongst(a," +
		            ",".join( [ i_roi
		                        for ( _, i_roi ) in conf[ "ana" ][ "vl_rois" ]
		                      ]
		                    ) +
		            ")'"
		          )

		vl_cmd = [ "3dcalc",
		           "-a", vl_path,
		           "-prefix", vl_sub_path,
		           "-expr", vl_expr,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( " ".join( vl_cmd ) )

		mask_path = paths.roi.mask_rois.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		# now we want to make a dataset from the non vis_loc rois
		dset_cmd = [ "ROI2dataset",
		             "-pad_to_node", "{nk:d}".format( nk = conf[ "subj" ][ "node_k" ][ hemi ] ),
		             "-prefix", mask_path,
		             "-overwrite",
		             "-input", " ".join( [ "{hemi:s}_{roi:s}.1D.roi".format( hemi = hemi, roi = roi )
		                                   for ( roi, _ ) in conf[ "ana" ][ "mask_rois" ]
		                                 ]
		                               )
		           ]

		fmri_tools.utils.run_cmd( " ".join( dset_cmd ) )

		# now we want to combine the vis_loc and the mask ROI datasets
		roi_path = paths.roi.rois.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		cmb_expr = "'a+(iszero(a)*b)'"

		cmb_cmd = [ "3dcalc",
		            "-a", vl_sub_path,
		            "-b", mask_path,
		            "-expr", cmb_expr,
		            "-overwrite",
		            "-prefix", roi_path
		          ]

		fmri_tools.utils.run_cmd( " ".join( cmb_cmd ) )

	os.chdir( start_dir )


def roi_xtr( conf, paths ):
	"""Extract PSC and statistics data from ROIs"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running ROI extraction..." )

	start_dir = os.getcwd()

	os.chdir( paths.roi.base.dir() )

	for hemi in [ "lh", "rh" ]:

		roi_path = paths.roi.rois.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		dset_path = paths.ana.psc.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		txt_path = paths.roi.psc.full( "_{hemi:s}.txt".format( hemi = hemi ) )

		mask_path = paths.ana.mask.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		# 3dmaskdump won't overwrite, so need to manually remove any previous file
		if os.path.exists( txt_path ):
			os.remove( txt_path )

		mask_expr = "'-a {m:s} -expr step(a)'".format( m = mask_path )

		xtr_cmd = [ "3dmaskdump",
		            "-mask", roi_path,
		            "-cmask", mask_expr,
		            "-noijk",
		            "-o", txt_path,
		            roi_path,
		            dset_path
		          ]

		fmri_tools.utils.run_cmd( " ".join( xtr_cmd ) )

	os.chdir( start_dir )


def task( conf, paths ):
	"""Analyses performance on the behavioural task"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running task analysis..." )


	# half-open interval
	run_time_bins = np.arange( start = 0,
	                           stop = conf[ "exp" ][ "run_full_len_s" ],
	                           step = conf[ "ana" ][ "task_bin_res_s" ]
	                         )[ :-1 ]

	n_bins_per_run = len( run_time_bins )

	# final dimension is ( task active, response registered, condition )
	task_data = np.zeros( ( conf[ "subj" ][ "n_runs" ],
	                        n_bins_per_run,
	                        3
	                      )
	                    )

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		task_path = paths.logs.task.full( "_{r:d}.npy".format( r = ( i_run + 1 ) ) )
		resp_path = paths.logs.resp.full( "_{r:d}.npy".format( r = ( i_run + 1 ) ) )
		seq_path = paths.logs.seq.full( "_{r:d}.npy".format( r = ( i_run + 1 ) ) )

		# ( event x ( time, digit, polarity, task ) )
		task = np.load( task_path )
		# ( response x ( letter, time ) )
		resp = np.load( resp_path )
		# we don't care about which key
		resp = np.array( [ resps[ "time" ] for resps in resp ] )
		# ( event x ( time, block, block type, coherence, orientation, contrast, generated ) )
		seq = np.load( seq_path )

		for ( i_bin, t_start ) in enumerate( run_time_bins ):

			# Q1: was there an active target during this bin?
			#  - find the active task event as the last event that has an onset time less than the start
			#    of the current bin
			i_task = np.where( t_start >= task[ :, 0 ] )[ 0 ][ -1 ]
			task_data[ i_run, i_bin, 0 ] = task[ i_task, 3 ]

			# Q2: did the subject respond during the bin?
			#  - find the first response prior to the current bin
			i_resp = np.where( t_start >= resp )[ 0 ]
			#  - only included if it is less than the bin distance away
			if ( i_resp.size > 0 ) and ( t_start - resp[ i_resp[ -1 ] ] ) < conf[ "ana" ][ "task_bin_res_s" ]:
				task_data[ i_run, i_bin, 1 ] = 1
			else:
				task_data[ i_run, i_bin, 1 ] = 0

			# Q3: what was the stimulus condition during this bin?
			i_seq = np.where( t_start >= seq[ :, 0 ] )[ 0 ]
			if ( i_seq.size > 0 ):
				task_data[ i_run, i_bin, 2 ] = seq[ i_seq[ -1 ], 2 ]

	# now it can be saved
	np.save( paths.task.data.full( ".npy" ), task_data )

	# ... before analysing performance
	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	task_perf = np.empty( ( conf[ "ana" ][ "task_perf_n_bins" ],
	                        n_cond
	                      )
	                    )
	task_perf.fill( np.NAN )

	for i_bin in xrange( conf[ "ana" ][ "task_perf_n_bins" ] ):

		bin_data = []

		for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

			run_data = task_data[ i_run, :, : ]

			# we lose some events because of the time shifting
			run_shft_data = np.empty( ( run_data.shape[ 0 ] - i_bin, run_data.shape[ 1 ] ) )
			run_shft_data.fill( np.NAN )

			if i_bin > 0:
				run_shft_data[ :, 0 ] = run_data[ :-i_bin, 0 ]
				run_shft_data[ :, 2 ] = run_data[ :-i_bin:, 2 ]

				run_shft_data[ :, 1 ] = run_data[ i_bin:, 1 ]

			else:
				run_shft_data = run_data

			bin_data.append( run_shft_data )

		bin_data = np.vstack( bin_data )

		task_perf[ i_bin, : ] = [ np.corrcoef( bin_data[ bin_data[ :, 2 ] == i_cond, :2 ].T )[ 1, 0 ]
		                          for i_cond in np.arange( 1, n_cond + 1 )
		                        ]

	assert( np.sum( np.isnan( task_perf ) ) == 0 )

	np.savetxt( paths.task.perf.full( ".txt" ), task_perf )

	os.chdir( paths.task.data.dir() )

	# now to save in NIML format
	# save each condition separately
	for i_cond in xrange( n_cond ):

		cond_data = task_perf[ :, i_cond ]

		cond_path = paths.task.data.full( "_{c:d}".format( c = i_cond ) )

		cond_sub_paths = cond_path + "_sub"

		fmri_tools.utils.array_to_niml( cond_data, cond_path, cond_sub_paths )
