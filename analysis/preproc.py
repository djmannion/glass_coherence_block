"""
Set of routines to pre-process the fMRI data for the Glass patterns block
design fMRI experiment.
"""

from __future__ import division

import os, os.path
import logging

import numpy as np

import fmri_tools.preproc, fmri_tools.utils

import glass_coherence_block.exp.run


def convert( conf, paths ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running conversion..." )

	# aggregate the dicom directories
	raw_paths = ( paths.func.raws +
	              [ paths.fmap.mag_raw, paths.fmap.ph_raw ]
	            )

	# aggregate the images paths
	img_paths = ( paths.func.origs +
	              [ paths.fmap.mag, paths.fmap.ph ]
	            )

	for ( raw_path, img_path ) in zip( raw_paths, img_paths ):

		fmri_tools.preproc.dcm_to_nii( raw_path.full(),
		                               img_path.full( ".nii" ),
		                               reorient_dim = conf[ "acq" ][ "ras" ]
		                             )


	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( [ img_path.full( ".nii" )
	                                             for img_path in img_paths
	                                           ]
	                                         )
	      )

	# files to go into the summary
	summ_paths = [ orig.full() for orig in paths.func.origs ]

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( summ_paths,
	                                      paths.summ.orig.full(),
	                                    )


def mot_correct( conf, paths ):
	"""Performs motion correction"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running motion correction..." )

	orig_paths = [ orig_path.full() for orig_path in paths.func.origs ]
	corr_paths = [ corr_path.full() for corr_path in paths.func.corrs ]

	i_base = conf[ "subj" ][ "mot_base" ] - 1
	base_path = "{fname:s}[0]".format( fname = paths.func.origs[ i_base ].full( ".nii" ) )

	fmri_tools.preproc.mot_correct( orig_paths = orig_paths,
	                                corr_paths = corr_paths,
	                                base_path = base_path,
	                                mc_path = paths.summ.motion.full( ".txt" )
	                              )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( corr_paths,
	                                      paths.summ.corr.full()
	                                    )


def fieldmaps( conf, paths ):
	"""Prepare the fieldmaps"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running fieldmap preparation..." )

	fmri_tools.preproc.make_fieldmap( mag_path = paths.fmap.mag.full(),
	                                  ph_path = paths.fmap.ph.full(),
	                                  fmap_path = paths.fmap.fmap.full(),
	                                  delta_te_ms = conf[ "acq" ][ "delta_te_ms" ]
	                                )

def unwarp( conf, paths ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running distortion correction..." )

	for ( corr_path, uw_path ) in zip( paths.func.corrs, paths.func.uws ):

		fmri_tools.preproc.unwarp( epi_path = corr_path.full(),
		                           fmap_path = paths.fmap.fmap.full(),
		                           uw_path = uw_path.full(),
		                           dwell_ms = conf[ "acq" ][ "dwell_ms" ],
		                           uw_direction = conf[ "acq" ][ "ph_encode_dir" ],
		                           pass_nocheck = False
		                         )

	uw_files = [ uw.full() for uw in paths.func.uws ]

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( uw_files, paths.summ.uw.full() )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( uw_files, paths.summ.mean.full() )


def sess_reg( conf, paths ):
	"""Coregisters the session anatomical with its mean functional"""

	logger = logging.getLogger( __name__ )

	logger.info( "Running registration..." )

	if conf[ "subj" ][ "extra_al_params" ]:
		extra_al_params = conf[ "subj" ][ "extra_al_params" ]
	else:
		extra_al_params = None

	fmri_tools.preproc.img_reg( reg_dir = paths.reg.base.full(),
	                            base_file = paths.reg.mean.file( "+orig" ),
	                            mov_file = paths.reg.anat_ref.file( "+orig" ),
	                            extra_al_params = extra_al_params
	                          )



def vol_to_surf( conf, paths ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	logger = logging.getLogger( __name__ )
	logger.info( "Running volume to surface projection..." )

	start_dir = os.getcwd()

	for ( uw_file, surf_file, run_dir ) in zip( paths.func.uws, paths.func.surfs, paths.func.runs ):

		os.chdir( run_dir.full() )

		for hemi in [ "lh", "rh" ]:

			surf_cmd = [ "3dVol2Surf",
			             "-spec", paths.reg.spec.full( "_{hemi:s}.spec".format( hemi = hemi ) ),
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "15",
			             "-f_index", "nodes",
			             "-sv", paths.reg.anat_reg.full( "+orig" ),
			             "-grid_parent", uw_file.full( ".nii" ),
			             "-out_niml", surf_file.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( " ".join( surf_cmd ) )

	os.chdir( start_dir )


def design_prep( conf, paths ):
	"""Prepares the designs for GLM analysis"""

	# load a dictionary that contains which index in the run sequence refers to
	# what aspect of the experiment
	seq_info = glass_coherence_block.exp.run.get_seq_ind()

	# number of conditions corresponds to the number of coherence levels (we dont
	# count fixation as a condition)
	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	# we want to write out a time file for each condition, where every row is a
	# run and the columns are the block onset times
	run_files = [ open( paths.ana.stim_times.full( "_{c:d}.txt".format( c = cond_num ) ),
	                    "w"
	                  )
	              for cond_num in np.arange( 1, n_cond + 1 )
	            ]

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		run_times = []
		run_conds = []

		# load the sequence info for this run, which is output by the experiment
		# script
		run_seq = np.load( paths.logs.seq.full( "_{run:d}.npy".format( run = ( i_run + 1 ) ) ) )

		( n_evt, _ ) = run_seq.shape

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
			run_files[ i_evt_file ].write( "{et:.5f}\t".format( et = evt_time ) )

		# finished this run, so add a newline to every condition's file
		_ = [ run_file.write( "\n" ) for run_file in run_files ]

	# finished, so close all the open files
	_ = [ run_file.close() for run_file in run_files ]
