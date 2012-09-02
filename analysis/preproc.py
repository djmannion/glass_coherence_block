"""
Set of routines to pre-process the fMRI data for the Glass patterns block
design fMRI experiment.
"""

from __future__ import division

import os.path
import tempfile

import numpy as np

import fmri_tools.preproc, fmri_tools.utils

import glass_coherence_block.exp.run


def convert( paths, conf ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# aggregate the dicom directories
	raw_dirs = [ paths[ "func" ][ "raw_dirs" ] +
	             paths[ "fmap" ][ "raw_mag_dirs" ] +
	             paths[ "fmap" ][ "raw_ph_dirs" ]
	           ]

	# aggregate the output directories
	nii_dirs = [ paths[ "func" ][ "run_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ]
	           ]

	# aggregate the images paths
	img_paths = [ paths[ "func" ][ "orig_files" ] +
	              paths[ "fmap" ][ "mag_files" ] +
	              paths[ "fmap" ][ "ph_files" ]
	            ]

	# pull out the filenames
	img_names = [ os.path.split( img_path )[ 1 ]
	              for img_path in img_paths
	            ]

	for i_dir in xrange( len( raw_dirs ) ):

		fmri_tools.preproc.dcm_to_nii( raw_dirs[ i_dir ],
		                               nii_dirs[ i_dir ],
		                               img_names[ i_dir ],
		                               reorient_dim = conf[ "acq" ][ "ras" ],
		                               log_path = paths[ "summ" ][ "log_file" ]
		                             )


	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_img_paths = [ "%s.nii" % img_path for img_path in img_paths ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_img_paths ) )

	# files to go into the summary
	summ_paths = paths[ "func" ][ "orig_files" ]

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( summ_paths,
	                                      paths[ "summ" ][ "orig_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
	                                    )


def run_masks( paths, conf ):
	"""Make a brainmask for each run and compute the inverse variance within."""

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		# the run's original file
		base_file = "%s.nii" % paths[ "func" ][ "orig_files" ][ i_run ]
		# output mask file
		mask_file = "%s.nii" % paths[ "func" ][ "mask_files" ][ i_run ]

		# the 'SI' flag gets rid of the cerebellum
		mask_cmd = [ "3dAutomask",
		             "-prefix", mask_file,
		             "-SI", "%.3f" % conf[ "subj" ][ "mask_z_mm" ],
		             "-overwrite",
		             base_file
		           ]

		fmri_tools.utils.run_cmd( mask_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# to compute the variance, we use '3dTstat' to first compute the standard
		# deviation. we don't care about this though, so we just save it as a temp
		# file so it is deleted after
		std_tmp = tempfile.NamedTemporaryFile( suffix = ".nii" )

		std_cmd = [ "3dTstat",
		            "-stdev",
		            "-overwrite",
		            "-prefix", std_tmp.name,
		            base_file
		          ]

		fmri_tools.utils.run_cmd( std_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		ivar_file = "%s.nii" % paths[ "func" ][ "ivar_files" ][ i_run ]

		# compute the inverse variance, with implicit masking
		# i got this trick from an AFNI forum post by Bob Cox
		ivar_cmd = [ "3dcalc",
		             "-a", std_tmp.name,
		             "-b", mask_file,
		             "-expr", "b/(a*a)",
		             "-overwrite",
		             "-prefix", ivar_file
		           ]

		fmri_tools.utils.run_cmd( ivar_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )


def motion_correct( paths, conf ):
	"""Perform motion correction"""

	# the index to the run that will be the 'base', to be corrected to
	# it is stored in one-based, hence the minus 1 to get to an index
	i_mc_base = conf[ "subj" ][ "mot_base" ] - 1

	mc_base = "%s.nii[0]" % paths[ "func" ][ "orig_files" ][ i_mc_base ]

	mc_params = []

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		orig_file = paths[ "func" ][ "orig_files" ][ i_run ]
		corr_file = paths[ "func" ][ "corr_files" ][ i_run ]
		ivar_file = paths[ "func" ][ "ivar_files" ][ i_run ]

		# because we want to aggregate the motion correction files over the whole
		# session, we only store the individual run correction temporarily
		mc_txt = tempfile.NamedTemporaryFile()

		mc_cmd = [ "3dvolreg",
		           "-twopass",
		           "-prefix", "%s.nii" % corr_file,
		           "-1Dfile", mc_txt.name,
		           "-overwrite",
		           "-weight", "%s.nii[0]" % ivar_file,
		           "-base", mc_base,
		           "-zpad", "5",
		           "-heptic",  # fourier can cause ringing artefacts
		           "%s.nii" % orig_file
		         ]

		fmri_tools.utils.run_cmd( mc_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# deal with the motion estimates
		mc_params.append( np.loadtxt( mc_txt.name ) )

	# concatenate the motion estimates over time
	mc_params = np.vstack( mc_params )

	np.savetxt( paths[ "summ" ][ "mot_est_file" ], mc_params )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( paths[ "func" ][ "corr_files" ],
	                                      paths[ "summ" ][ "corr_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
	                                    )


def fieldmaps( paths, conf ):
	"""Prepare the fieldmaps"""

	for i_fmap in xrange( conf[ "subj" ][ "n_fmaps" ] ):

		fmri_tools.preproc.make_fieldmap( paths[ "fmap" ][ "mag_files" ][ i_fmap ],
		                                  paths[ "fmap" ][ "ph_files" ][ i_fmap ],
		                                  paths[ "fmap" ][ "fmap_files" ][ i_fmap ],
		                                  conf[ "acq" ][ "delta_te_ms" ],
		                                  log_path = paths[ "summ" ][ "log_file" ]
		                                )


def unwarp( paths, conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	func_fmap = paths[ "fmap" ][ "fmap_files" ][ 0 ]

	# motion-corrected images (input)
	func_corr = paths[ "func" ][ "corr_files" ]
	# unwarped images (output)
	func_uw = paths[ "func" ][ "uw_files" ]

	for i_run in xrange( len( func_corr )  ):

		fmri_tools.preproc.unwarp( func_corr[ i_run ],
		                           func_fmap,
		                           func_uw[ i_run ],
		                           conf[ "acq" ][ "dwell_ms" ],
		                           conf[ "acq" ][ "ph_encode_dir" ],
		                           log_path = paths[ "summ" ][ "log_file" ]
		                         )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( func_uw,
	                               paths[ "summ" ][ "mean_file" ],
	                               log_path = paths[ "summ" ][ "log_file" ]
	                             )

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( func_uw,
	                                      paths[ "summ" ][ "uw_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
	                                    )


def surf_reg( paths, conf ):
	"""Coregisters an anatomical with the SUMA reference"""

	base_anat = paths[ "reg" ][ "anat" ]
	reg_anat = paths[ "reg" ][ "reg_anat" ]

	# use the mean to represent the functional images
	base_func = "%s.nii" % paths[ "summ" ][ "mean_file" ]

	coreg_cmd = [ "3dAllineate",
	              "-base", base_func,
	              "-source", base_anat,
	              "-prefix", reg_anat,
	              "-cost", "nmi",  # normalised mutual info cost function
	              "-master", "SOURCE",
	              "-warp", "shift_rotate",  # only rigid transforms
	              "-onepass",  # false minima if it is allowed to wander
	              "-verb"
	            ]

	# pass the algorithm three translation parameters to get it close to the
	# reference anatomy
	for ( i_nudge, nudge_val ) in enumerate( conf[ "subj" ][ "nudge_vals" ] ):

		coreg_cmd.extend( [ "-parini",
		                    "%d" % ( i_nudge + 1 ),
		                    "%.3f" % nudge_val
		                  ]
		                )

	fmri_tools.utils.run_cmd( coreg_cmd,
	                          env = fmri_tools.utils.get_env(),
	                          log_path = paths[ "summ" ][ "log_file" ]
	                        )


def vol_to_surf( paths, conf ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	# this puts some output in the working directory, so change to where we want
	# to save stuff
	start_dir = os.getcwd()

	# images to project (unwarped)
	vol_files = paths[ "func" ][ "uw_files" ]
	# surface files to write
	surf_files = paths[ "func" ][ "surf_files" ]

	# number of increments to divide the surface interval
	f_steps = 15
	# how to deal with multiple nodes lying on a given voxel
	f_index = "nodes"

	for ( vol_file, surf_file ) in zip( vol_files, surf_files ):

		( file_dir, _ ) = os.path.split( vol_file )

		os.chdir( file_dir )

		for hemi in [ "lh", "rh" ]:

			surf_cmd = [ "3dVol2Surf",
			             "-spec", "%s%s.spec" % ( paths[ "reg" ][ "spec" ], hemi ),
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "%d" % f_steps,
			             "-f_index", f_index,
			             "-sv", paths[ "reg" ][ "reg_anat" ],
			             "-grid_parent", "%s.nii" % vol_file,
			             "-out_niml", "%s_%s.niml.dset" % ( surf_file, hemi ),
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( surf_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	os.chdir( start_dir )


def design_prep( paths, conf ):
	"""Prepares the designs for GLM analysis"""

	# load a dictionary that contains which index in the run sequence refers to
	# what aspect of the experiment
	seq_info = glass_coherence_block.exp.run.get_seq_ind()

	# number of conditions corresponds to the number of coherence levels (we dont
	# count fixation as a condition)
	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	# we want to write out a time file for each condition, where every row is a
	# run and the columns are the block onset times
	run_files = [ open( "%s%d.txt" % ( paths[ "ana" ][ "time_files" ],
	                                   cond_num
	                                 ),
	                    "w"
	                  )
	              for cond_num in np.arange( 1, n_cond + 1 )
	            ]

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		run_times = []
		run_conds = []

		# load the sequence info for this run, which is output by the experiment
		# script
		run_seq = np.load( "%s%d.npy" % ( paths[ "log" ][ "seq_base" ], i_run + 1 ) )

		( n_evt, n_params ) = run_seq.shape

		# loop through each event in the sequence
		for i_evt in xrange( n_evt ):

			# pull out which block this and the previous event comes from
			curr_block_num = run_seq[ i_evt, seq_info[ "block_num" ] ]
			prev_block_num = run_seq[ i_evt - 1, seq_info[ "block_num" ] ]

			# it is a 'transition' event if it has a different block number
			is_transition = curr_block_num != prev_block_num

			# we only care about the transition blocks
			if is_transition:

				# pull out the onset time
				start_time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

				# and the condition number, cast to an int so it can be used as an
				# index
				cond = int( run_seq[ i_evt, seq_info[ "block_type" ] ] )

				# store the start time and condition type in the relevent list
				run_times.append( start_time_s )
				run_conds.append( cond )

		# finished iterating through the events, now convert from lists to arrays
		run_times = np.array( run_times )
		run_conds = np.array( run_conds )

		# test for valid events as those that lie within the run duration
		valid_evts = np.logical_and( run_times >= 0,
		                             run_times < conf[ "exp" ][ "run_len_s" ]
		                           )

		# and don't include fixation
		valid_evts = np.logical_and( valid_evts, run_conds > 0 )

		# restrict our data to the valid events
		run_times = run_times[ valid_evts ]
		run_conds = run_conds[ valid_evts ]

		for ( evt_time, evt_cond ) in zip( run_times, run_conds ):

			# the conditions are one-based
			i_evt_file = evt_cond - 1

			# write the onset time
			run_files[ i_evt_file ].write( "%.5f\t" % evt_time )

		# finished this run, so add a newline to every condition's file
		_ = [ run_file.write( "\n" ) for run_file in run_files ]

	# finished, so close all the open files
	_ = [ run_file.close() for run_file in run_files ]


	# POLYNOMIALS
	# ---

	# total number of volumes
	n_tot_vol = conf[ "exp" ][ "run_full_len_s" ] / conf[ "acq" ][ "tr_s" ]

	# number of volumes to reject at the start and end
	n_pre_vol = conf[ "exp" ][ "pre_len_s" ] / conf[ "acq" ][ "tr_s" ]
	n_post_vol = conf[ "exp" ][ "post_len_s" ] / conf[ "acq" ][ "tr_s" ]

	# number of volumes left after rejection
	n_valid_vol = n_tot_vol - n_pre_vol - n_post_vol

	# compute the polynomial timecourses
	run_trends = fmri_tools.utils.legendre_poly( conf[ "ana" ][ "poly_ord" ],
	                                             int( n_valid_vol ),
	                                             pre_n = int( n_pre_vol ),
	                                             post_n = int( n_post_vol )
	                                           )

	assert( run_trends.shape[ 0 ] == n_tot_vol )

	# need to have a set of trends for each run, zeroed elsewhere
	bl_trends = np.zeros( ( n_tot_vol  * conf[ "subj" ][ "n_runs" ],
	                        conf[ "ana" ][ "poly_ord" ] * conf[ "subj" ][ "n_runs" ]
	                      )
	                    )

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		i_row_start = i_run * run_trends.shape[ 0 ]
		i_row_end = i_row_start + run_trends.shape[ 0 ]

		i_col_start = i_run * run_trends.shape[ 1 ]
		i_col_end = i_col_start + run_trends.shape[ 1 ]

		bl_trends[ i_row_start:i_row_end, i_col_start:i_col_end ] = run_trends

	np.savetxt( paths[ "ana" ][ "bl_poly" ], bl_trends )
