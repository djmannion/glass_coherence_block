"""
Set of routines to pre-process the fMRI data for the Glass coherence block
design fMRI experiment.
"""

from __future__ import division

import os.path

import nipy
import numpy as np
import scipy.io

import glass_coherence_block.exp.run
import fmri_tools.preproc, fmri_tools.utils

def convert( paths ):
	"""Converts the functionals and fieldmaps from dicom to nifti.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence_block.config'.

	"""

	# aggregate the dicom directories
	dcm_dirs = ( paths[ "func" ][ "dcm_dirs" ] +
	             paths[ "fmap" ][ "dcm_mag_dirs" ] +
	             paths[ "fmap" ][ "dcm_ph_dirs" ]
	           )

	# aggregate the output directories
	nii_dirs = ( paths[ "func" ][ "dirs" ] +
	             paths[ "fmap" ][ "dirs" ] +
	             paths[ "fmap" ][ "dirs" ]
	           )

	# aggregate the images paths
	img_paths = ( paths[ "func" ][ "orig" ] +
	              paths[ "fmap" ][ "mag" ] +
	              paths[ "fmap" ][ "ph" ]
	            )

	# pull out the filenames
	img_names = [ os.path.split( img_path )[ 1 ]
	              for img_path in img_paths
	            ]

	# do the DCM -> NII conversion
	map( fmri_tools.preproc.dcm_to_nii,
	     dcm_dirs,
	     nii_dirs,
	     img_names
	   )

	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_img_paths = [ "%s.nii" % img_path
	                   for img_path in img_paths
	                 ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_img_paths ) )

	# make a summary image from the functional images
	fmri_tools.preproc.gen_sess_summ_img( paths[ "func" ][ "orig" ],
	                                      paths[ "func" ][ "orig_summ" ]
	                                    )


def st_motion_correct( paths, conf, subj_conf ):
	"""Performs slice-timing and motion correction.
	"""

	# get the order of runs to pass to the correction algorithm
	# this is done because the algorithm realigns all to the first entry, which
	# normally corresponds to the first run acquired, but we might want them to
	# be aligned with a different run - one closer to a fieldmap, for example
	run_order = subj_conf[ "run_st_mot_order" ]

	# reorder the paths
	# (the -1 is because the runs are specified in subj_conf in a one-based
	# index; ie. run 1 is the first run)
	orig_paths = [ paths[ "func" ][ "orig" ][ i_run - 1 ]
	               for i_run in run_order
	             ]
	corr_paths = [ paths[ "func" ][ "corr" ][ i_run - 1 ]
	               for i_run in run_order
	             ]

	# pull out the important information from the config
	slice_order = conf[ "acq" ][ "slice_order" ]
	tr_s = conf[ "acq" ][ "tr_s" ]
	slice_info = ( conf[ "acq" ][ "slice_axis" ],
	               conf[ "acq" ][ "slice_acq_dir" ]
	             )

	# run the motion correction algorithm (slow)
	motion_est = fmri_tools.preproc.correct_st_motion( orig_paths,
	                                                   corr_paths,
	                                                   slice_order,
	                                                   tr_s,
	                                                   slice_info,
	                                                 )

	# save the estimated motion parameters
	np.save( paths[ "func" ][ "motion_estimates" ],
	         arr = motion_est
	       )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( corr_paths,
	                                      paths[ "func" ][ "corr_summ" ]
	                                    )


def fieldmaps( paths, conf, subj_conf ):
	"""Prepare the fieldmaps.
	"""

	# duplicate the delta TE for each fieldmap acquired
	delta_te_ms = ( [ conf[ "acq" ][ "delta_te_ms" ] ] *
	                subj_conf[ "n_fmaps" ]
	              )

	map( fmri_tools.preproc.make_fieldmap,
	     paths[ "fmap" ][ "mag" ],
	     paths[ "fmap" ][ "ph" ],
	     paths[ "fmap" ][ "fmap" ],
	     delta_te_ms
	   )


def unwarp( paths, conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	# combine the experiment and localiser functional info
	func_corr = paths[ "func" ][ "corr" ]
	func_fmap = paths[ "func" ][ "fmap" ]
	func_uw = paths[ "func" ][ "uw" ]

	# duplicate the dwell time and phase encode direction for each image
	dwell_ms = [ conf[ "acq" ][ "dwell_ms" ] ] * len( func_corr )
	ph_encode_dir = [ conf[ "acq" ][ "ph_encode_dir" ] ] * len( func_corr )

	# perform the unwarping
	map( fmri_tools.preproc.unwarp,
	     func_corr,
	     func_fmap,
	     func_uw,
	     dwell_ms,
	     ph_encode_dir
	   )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( func_uw,
	                               paths[ "func" ][ "mean" ]
	                             )

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( func_uw,
	                                      paths[ "func" ][ "uw_summ" ]
	                                    )


def make_roi_images( paths, conf ):
	"""Converts the ROI matlab files to nifti images in register with the
	subject's anatomical.
	"""

	n_rois = len( conf[ "ana" ][ "rois" ] )

	# python nifti loader needs the extension
	anat_path = [ "%s.nii" % paths[ "anat" ][ "anat" ] ] * n_rois

	roi_paths = [ "%s.nii" % roi_path
	              for roi_path in paths[ "roi" ][ "orig" ]
	            ]

	roi_ax = [ conf[ "ana" ][ "roi_ax" ] ] * n_rois
	roi_ax_order = [ conf[ "ana" ][ "roi_ax_order" ] ] * n_rois

	# convert the ROI mat files to nifti images
	map( fmri_tools.preproc.roi_to_nii,
	     paths[ "roi" ][ "mat" ],
	     anat_path,
	     roi_paths,
	     roi_ax,
	     roi_ax_order
	   )


def prepare_rois( paths, conf ):
	"""Extracts and writes the coordinates of each ROI.
	"""

	for ( i_roi, roi_name ) in enumerate( conf[ "ana" ][ "rois" ] ):

		# load the ROI image
		roi = nipy.load_image( "%s.nii" % paths[ "roi" ][ "rs" ][ i_roi ]
		                     ).get_data()

		# get rid of any NaNs
		roi[ np.isnan( roi) ] = 0

		# extract the coordinate list
		coords = np.nonzero( roi )

		# and save
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "coords" ], roi_name ),
		         arr = coords
		       )


def form_vtcs( paths, conf, subj_conf ):
	"""Extracts the voxel time courses for each voxel in each ROI
	"""

	n_full_vols_per_run = ( conf[ "exp" ][ "run_full_len_s" ] /
	                        conf[ "acq" ][ "tr_s" ]
	                      )

	n_vols_per_run = ( conf[ "exp" ][ "run_len_s" ] /
	                   conf[ "acq" ][ "tr_s" ]
	                 )

	n_start_cull_vols = ( conf[ "exp" ][ "pre_len_s" ] /
	                      conf[ "acq" ][ "tr_s" ]
	                    )

	n_end_cull_vols = ( conf[ "exp" ][ "post_len_s" ] /
	                    conf[ "acq" ][ "tr_s" ]
	                  )

	valid_vol_range = np.arange( n_start_cull_vols,
	                             n_full_vols_per_run - n_end_cull_vols
	                           ).astype( "int" )

	assert( len( valid_vol_range ) == n_vols_per_run )


	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the coordinates for the roi
		coords = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "coords" ], roi_name ) )

		n_voxels = coords.shape[ 1 ]

		# initialise the vtc
		vtc = np.empty( ( n_full_vols_per_run,
		                  subj_conf[ "n_runs" ],
		                  n_voxels
		                )
		              )

		# fill with NaNs, to be safe
		vtc.fill( np.NAN )

		# loop through each unwarped image (run) file
		for ( i_run, run_path ) in enumerate( paths[ "func" ][ "uw" ] ):

			# load the run file
			run_img = nipy.load_image( "%s.nii" % run_path ).get_data()

			# iterate through each voxel in the roi
			for i_voxel in xrange( n_voxels ):

				# extract the voxel data (timecourse) at the voxel coordinate
				vox_data = run_img[ coords[ 0, i_voxel ],
				                    coords[ 1, i_voxel ],
				                    coords[ 2, i_voxel ],
				                    :
				                  ]

				# store the voxel data
				vtc[ :, i_run, i_voxel ] = vox_data

		# discard the unwanted volumes
		vtc = vtc[ valid_vol_range, :, : ]

		# check that it has been filled up correctly
		assert( np.sum( np.isnan( vtc ) ) == 0 )

		# save the vtc
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "vtc" ], roi_name ),
		         arr = vtc
		       )


def roi_vtc_cull( paths, conf ):
	"""Culls voxels from each ROI based on their mean-normalised variance."""

	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the coordinates
		coords = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "coords"], roi_name ) )

		# load the vtc
		vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "vtc" ], roi_name ) )

		# do a quick check to make sure the coords have the expected shape
		assert( coords.shape[ 0 ] == 3 )
		assert( coords.shape[ 1 ] == vtc.shape[ -1 ] )

		# this is a logical matrix showing if the voxel passes the selection
		# criteria for each run
		retain = np.empty( ( vtc.shape[ 1:3 ] ) )

		cull_p = conf[ "ana" ][ "cull_prop" ]

		for i_run in xrange( vtc.shape[ 1 ] ):

			run_vtc = vtc[ :, i_run, : ]

			retain[ i_run, : ] = fmri_tools.preproc.roi_vox_trim_var( run_vtc,
			                                                          cull_p
			                                                        )

		# to be retained, has to be retained for all runs
		vox_to_retain = np.all( retain, axis = 0 )

		# do the cullin'
		vtc = vtc[ :, :, vox_to_retain ]
		coords = coords[ :, vox_to_retain ]

		# and re-save
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "coords_sel" ], roi_name ),
		         coords
		       )

		# save the vtc
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_sel" ], roi_name ),
		         vtc
		       )


def get_evt_info( paths, conf, subj_conf ):
	"""Extracts and writes the stimulus event information from the session log.
	"""

	# each event will be appended to this
	evt_info = []

	# so we know what each indice in the log file means
	i_evt_log = glass_coherence_block.exp.run.get_seq_ind()

	for i_run in xrange( subj_conf[ "n_runs" ] ):

		# load the run sequence
		log = np.load( os.path.join( paths[ "design" ][ "log_dir" ],
		                             "%s_glass_coherence_block_seq_%d.npy" % (
		                             subj_conf[ "subj_id" ], i_run + 1 )
		                           )
		             )

		# iterate through each event in the log
		for i_evt in xrange( log.shape[ 0 ] ):

			evt = log[ i_evt, : ]

			# pull out the event's onset time
			evt_on_s = evt[ i_evt_log[ "time_s" ] ]
			# adjust the onset time to take into account the time at the beginning
			# that we cull
			evt_on_s -= conf[ "exp" ][ "pre_len_s" ]

			evt_off_s = evt_on_s + conf[ "exp" ][ "evt_stim_s" ]

			evt_cond = evt[ i_evt_log[ "block_type" ] ]

			evt_ori = np.round( evt[ i_evt_log[ "ori" ] ] )

			# don't include it if it is a blank condition
			if evt_cond > 0:

				if evt_cond == 1:
					if evt_ori==0:
						evt_c = 1
					elif evt_ori==90:
						evt_c = 2
				elif evt_cond == 2:
					if evt_ori==0:
						evt_c = 3
					elif evt_ori==90:
						evt_c = 4
				elif evt_cond == 3:
					if evt_ori==0:
						evt_c = 5
					elif evt_ori==90:
						evt_c = 6
				if evt_cond == 4:
					if evt_ori==0:
						evt_c = 7
					elif evt_ori==90:
						evt_c = 8

				evt_info.append( [ ( i_run + 1 ),
				                   evt_on_s,
				                   evt_off_s,
				                   evt_c
				                 ]
				               )

	evt_info = np.array( evt_info )

	np.save( paths[ "design" ][ "evt_info" ], evt_info )


def avg_vtcs( paths, conf ):
	"""Averages the timecourses over all the (selected) voxels in a ROI"""

	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the (post voxel selection) vtc
		vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_sel" ],
		                               roi_name
		                             )
		             )

		# average over voxels
		vtc = np.mean( vtc, axis = 2 )

		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_avg" ], roi_name ),
		         vtc
		       )
