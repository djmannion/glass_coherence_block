"""Configuration for the Glass coherence block design fMRI experiment.
"""

from __future__ import division

import os
import string

import numpy as np

import fmri_tools.utils


def get_conf():
	"""Overall experiment configuration.

	Returns
	-------
	conf : dict, with items:
		exp : holds overall experiment configurations
		stim : stimulus configuration
		task : behavioural task configuration
		acq : acquisition configuration

	"""

	conf = { "exp" : _get_exp_conf(),
	         "stim" : _get_stim_conf(),
	         "task" : _get_task_conf(),
	         "acq" : _get_acq_conf(),
	         "ana" : _get_ana_conf()
	       }

	return conf


def _get_ana_conf():
	"""
	"""

	ana_conf = { "rois" : ( "V1", "V2", "V3", "V3AB", "hV4" ),
	             "roi_ax" : ( 2, 1, 0 ),
	             "roi_ax_order" : ( 1, -1, -1 ),
	             "poly_ord" : 4
	           }

	return ana_conf


def _get_stim_conf():
	"""Get the stimulus configuration.

	Specifies the stimulus configuration and parameters for the Glass coherence
	block design fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------------

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	stim_conf = {}

	stim_conf[ "diam_deg" ] = 14.4

	# gives a density of 25 dots / deg ^2
	stim_conf[ "n_dipoles" ] = 2592

	# this is 6*sigma
	stim_conf[ "dot_size_deg" ] = 0.15

	stim_conf[ "pole_sep_deg" ] = 0.14

	stim_conf[ "ori_deg" ] = np.array( ( 0, 90 ) )

	stim_conf[ "coh_levels" ] = [ 0, 0.33, 0.66, 1 ]

	stim_conf[ "mask_in_deg" ] = 1.5
	stim_conf[ "mask_out_deg" ] = 7.2
	stim_conf[ "mask_fringe_deg" ] = 0.75

	return stim_conf


def _get_exp_conf():
	"""Gets the experiment configuration.

	Specifies the experiment configuration and parameters for the Glass coherence
	block design fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------
	n_blocks : scalar integer
		Number of blocks per run.
	block_len_s : scalar float
		Length of each block, in seconds.
	run_dur_s : scalar float
		Length of each run, in seconds.
	n_evt_per_block : scalar integer
		Number of events per block.
	n_evt_per_run : scalar integer
		Number of events per run.
	evt_len_s : scalar float
		Length of each 'event' within a block, in seconds.
	evt_stim_s : scalar float
		Length of stimulus presentation within each event.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	exp_conf = {}

	exp_conf[ "n_blocks" ] = 21

	exp_conf[ "block_len_s" ] = 16.0

	exp_conf[ "pre_len_s" ] = 6
	exp_conf[ "post_len_s" ] = 0

	exp_conf[ "run_len_s" ] = exp_conf[ "n_blocks" ] * exp_conf[ "block_len_s" ]

	exp_conf[ "run_full_len_s" ] = ( exp_conf[ "run_len_s" ] +
	                                 exp_conf[ "pre_len_s" ] +
	                                 exp_conf[ "post_len_s" ]
	                               )


	exp_conf[ "evt_len_s" ] = 1

	exp_conf[ "n_evt_per_run" ] = int( exp_conf[ "run_len_s" ] /
	                                   exp_conf[ "evt_len_s" ]
	                                 )

	exp_conf[ "evt_stim_s" ] = 0.75

	return exp_conf


def _get_task_conf():
	"""Gets the task configuration.

	Specifies the configuration for the behavioural task in the Glass coherence
	block design fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	task_conf = {}

	task_conf[ "set" ] = np.arange( 10 )

	task_conf[ "polarity" ] = ( -1, +1 )

	task_conf[ "rate_hz" ] = 3.0

	return task_conf


def _get_acq_conf():
	"""Get the acquisition configuration.

	Specifies the acquisition configuration and parameters for the
	Glass coherence block design fMRI experiment.

	Returns
	-------
	monitor_name : string
		monitor configuration name.
	tr_s : float
		time-to-repetition (TR), in seconds.
	delta_te_ms : float
		echo time differences for the fieldmaps, in milliseconds.
	dwell_ms : float
		dwell time, in milliseconds.
	slice_order : array of int
		slice acqusition indices, where 0 is the first slice.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	acq_conf = {}

	acq_conf[ "monitor_name" ] = "UMN_7T"

	acq_conf[ "tr_s" ] = 2.0

	acq_conf[ "delta_te_ms" ] = 1.02

	acq_conf[ "dwell_ms" ] = 0.325

	acq_conf[ "slice_order" ] = fmri_tools.utils.get_slice_order( 36 )

	acq_conf[ "slice_axis" ] = 2

	acq_conf[ "slice_acq_dir" ] = 1

	acq_conf[ "ph_encode_dir" ] = "y-"

	return acq_conf


def get_subj_conf():
	"""Gets the configuration info for each subject.

	Returns
	-------
	subj_id : string
		Subject ID, in the Olman lab system.
	acq_date : string
		Acquisition date, in YYYYMMDD format.
	comments : string
		Any comments about the scanning session.

	Notes
	-----
	* Return values are within a dictionary, which is itself within a dictionary
	  indexed by subject ID.

	"""

	s1021 = { "subj_id" : "s1021",
	          "acq_date" : "20120522",
	          "comments" : "",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 )
	        }

	subj_conf = { "s1021" : s1021 }

	return subj_conf


def get_exp_paths():
	"""
	"""

	# hard code
	BASE_PATH = "/home/dmannion/GlassCoherenceBlock/"

	paths = {}

	paths[ "dir" ] = BASE_PATH

	paths[ "grp_dir" ] = os.path.join( paths[ "dir" ],
	                                   "group_data"
	                                 )

	paths[ "grp_task" ] = os.path.join( paths[ "grp_dir" ],
	                                    "grp_task.mat"
	                                  )

	return paths


def get_subj_paths( subj_id ):
	"""Gets the filesystem path and file structure for the experiment data

	Parameters
	----------
	subj_id : string
		Subject identification number, eg. s1021

	Returns
	-------
	A dictionary with the following path info
	[ "exp" ]
		Paths at the experiment level
	[ "subj" ]
		Paths at the subject level
	[ "dir" ]
		Base directory for the subject.
	[ "func" ]
		Paths for the functional data.
	[ "fmap" ]
		Paths for the fieldmap data.
	[ "anat" ]
		Paths for the anatomical data.
	[ "roi" ]
		Paths for the ROI data.
	[ "analysis" ]
		Paths for the data analysis.
	[ "task" ]
		Paths for the task data.
	[ "design" ]
		Path to the design matrix.
	"""

	def get_task_paths( subj_conf, subj_dir ):
		""" Returns the paths for the fixation task

		Parameters
		----------
		subj_conf: dict
			Subject configuration directory
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		task : dict of strings
			Dictionary with the following path info:
			[ "dir" ] :
				Base directory for the task.
			[ "task_info" ] :
				Base path for the task information data.
			[ "task_resp" ] :
				Base path for the task response data.
		"""

		task = {}

		task[ "dir" ] = os.path.join( subj_dir, "log" )


		return task


	def get_analysis_paths( subj_dir ):
		""" Returns the paths for the analysis files

		Parameters
		----------
		subj_conf: dict
			Subject configuration directory
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		analysis : dict of strings
			Dictionary with the following path info:
			[ "dir" ] :
				Base directory for the analysis.
			[ "vtc" ] :
				Base path for the voxel timecourses of each ROI.
			[ "coords" ] :
				Base path for the coordinates of each ROI.
		"""

		analysis = {}

		analysis[ "dir" ] = os.path.join( subj_dir, "analysis" )

		analysis[ "vtc" ] = os.path.join( analysis[ "dir" ], "vtc" )

		analysis[ "vtc_sel" ] = os.path.join( analysis[ "dir" ], "vtc_sel" )
		analysis[ "vtc_avg" ] = os.path.join( analysis[ "dir" ], "vtc_avg" )

		analysis[ "coords" ] = os.path.join( analysis[ "dir" ], "coords" )
		analysis[ "coords_sel" ] = os.path.join( analysis[ "dir" ], "coords_sel" )

		analysis[ "block" ] = os.path.join( analysis[ "dir" ], "block" )

		analysis[ "cond_diff" ] = os.path.join( analysis[ "dir" ], "cond_diff" )

		return analysis


	def get_func_paths( subj_conf, subj_dir ):
		""" Returns the paths for the functional acquisitions

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		func : dict of strings
			Dictionary with the following path info:
			[ "dir " ] :
				Base directory for the functional data.
			[ "dirs" ] :
				Directory of each functional run.
			[ "dcm_dirs" ] :
				Directory of each functional run's raw DICOM directory.
			[ "orig" ], [ "corr" ], [ "uw" ]
				Filenames for each functional run's original, corrected, and unwarped
				NII files.
			[ "fmap" ] :
				Filename of each functional run's fieldmap.
			[ "mean" ] :
				Filename of an image containing the mean of all functional data.
			[ "motion_estimates" ] :
				Filename of an array containing the estimated motion parameters.

		"""

		func = {}

		# base functional directory
		func[ "dir" ] = os.path.join( subj_dir, "func" )

		# directory of each functional run
		func[ "dirs" ] = [ os.path.join( func[ "dir" ],
		                                 "run%02d" % ( i_run + 1 )
		                               )
		                   for i_run in xrange( subj_conf[ "n_runs" ] )
		                 ]

		# directory of each functional run's raw dicoms
		func[ "dcm_dirs" ] = [ os.path.join( func_dir,
		                                     "raw"
		                                   )
		                       for func_dir in func[ "dirs" ]
		                     ]

		# path to each functional run's 4D nifti file
		func[ "orig" ] = [ os.path.join( func_dir,
		                                 "%s_glass_coherence_block_run_%02d-orig" %
		                                 ( subj_conf[ "subj_id" ],
		                                   i_run + 1
		                                 )
		                               )
		                   for ( i_run, func_dir ) in enumerate( func[ "dirs" ] )
		                 ]

		# path to each functional run's 4D nifti file, slice-time and motion
		# corrected
		func[ "corr" ] = [ ( "%s-corr" % orig )
		                   for orig in func[ "orig" ]
		                 ]

		# ... and unwarped (distortion corrected)
		func[ "uw" ] = [ ( "%s-uw" % corr )
		                 for corr in func[ "corr" ]
		               ]

		# path to the fieldmap for each functional run
		func[ "fmap" ] = [ string.replace( orig, "-orig", "-fmap" )
		                   for orig in func[ "orig" ]
		                 ]

		# path to a mean image of all functional images
		func[ "mean" ] = os.path.join( func[ "dir" ],
		                               "%s_glass_coherence_block-mean" %
		                               subj_conf[ "subj_id" ]
		                             )

		# path to a summary image of all functional images, after motion correction
		func[ "orig_summ" ] = os.path.join( func[ "dir" ],
		                                    "%s_glass_coherence_block-orig-summ" %
		                                    subj_conf[ "subj_id" ]
		                                  )

		# path to a mean image of all functional images, after motion correction
		func[ "corr_summ" ] = os.path.join( func[ "dir" ],
		                                    "%s_glass_coherence_block-corr-summ" %
		                                    subj_conf[ "subj_id" ]
		                                  )

		# path to a mean image of all functional images, after unwarping
		func[ "uw_summ" ] = os.path.join( func[ "dir" ],
		                                  "%s_glass_coherence_block-uw-summ" %
		                                  subj_conf[ "subj_id" ]
		                                )

		# path to numpy file containing the motion estimates from the correction
		# procedure
		func[ "motion_estimates" ] = os.path.join( func[ "dir" ],
		                                           "motion_estimates.npy"
		                                         )

		return func


	def get_fmap_paths( subj_conf, subj_dir ):
		""" Returns the paths for the fieldmap acquisitions

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		fmap : dict of paths
			[ "dir " ] :
				Base directory for the fieldmap data.
			[ "dirs" ] :
				Directory of each fieldmap acquisition.
			[ "dcm_mag_dirs" ], [ "dcm_ph_dirs" ]:
				Directory of each fieldmap's magnitude and phase raw DICOM directories.
			[ "mag" ], [ "ph" ] :
				Filenames for each fieldmap's magnitude and phase NII files.
			[ "fmap" ] :
				Filename of each fieldmap's final fieldmap.

		"""

		fmap = {}

		# base directory for the fieldmaps
		fmap[ "dir" ] = os.path.join( subj_dir, "fmap" )

		# directories of each of the fieldmaps acquired in the session
		fmap[ "dirs" ] = [ os.path.join( fmap[ "dir" ],
		                                 "f%d" % ( i_fmap + 1 )
		                               )
		                   for i_fmap in xrange( subj_conf[ "n_fmaps" ] )
		                 ]

		# directories holding the raw magnitude images for each fieldmap
		fmap[ "dcm_mag_dirs" ] = [ os.path.join( fmap_dir,
		                                         "mag-raw"
		                                       )
		                           for fmap_dir in fmap[ "dirs" ]
		                         ]

		# ... raw phase image directories
		fmap[ "dcm_ph_dirs" ] = [ os.path.join( fmap_dir,
		                                        "ph-raw"
		                                      )
		                          for fmap_dir in fmap[ "dirs" ]
		                        ]

		# paths to the nifti file of each fieldmap's magnitude image
		fmap[ "mag" ] = [ os.path.join( fmap_dir,
		                                "%s_glass_coherence_block_fmap_%d-mag" %
		                                ( subj_conf[ "subj_id" ], i_fmap + 1 )
		                              )
		                  for ( i_fmap, fmap_dir ) in enumerate( fmap[ "dirs" ] )
		                ]

		# ... phase image
		fmap[ "ph" ] = [ os.path.join( fmap_dir,
		                               "%s_glass_coherence_block_fmap_%d-ph" %
		                               ( subj_conf[ "subj_id" ], i_fmap + 1 )
		                             )
		                 for ( i_fmap, fmap_dir ) in enumerate( fmap[ "dirs" ] )
		               ]

		# paths to the nifti file of each fieldmap
		fmap[ "fmap" ] = [ os.path.join( fmap_dir,
		                                 "%s_glass_coherence_block_fmap_%d-fmap" %
		                                 ( subj_conf[ "subj_id" ], i_fmap + 1 )
		                               )
		                   for ( i_fmap, fmap_dir ) in enumerate( fmap[ "dirs" ] )
		                 ]

		return fmap


	def get_design_paths( subj_conf, subj_dir ):
		""" Returns the paths for the session design.
		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		design : dict of paths
			[ "matrix" ] :
				Path to the design matrix.
			[ "run_seq" ] :
				Paths to the sequence info for each run.

		"""

		design = {}

		design[ "log_dir" ] = os.path.join( subj_dir, "log" )

		design[ "design" ] = os.path.join( subj_dir,
		                                   "design.npy"
		                                 )

		return design


	def get_anat_paths( subj_conf, subj_dir ):
		""" Returns the paths for the anatomy

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		anat : dict of paths
			[ "dir" ] :
				Base directory for the anatomical data.
			[ "anat" ] :
				Filename of the anatomical image.

		"""

		anat = {}

		# base anatomical directory
		anat[ "dir" ] = os.path.join( subj_dir, "anat" )

		# path to the anatomical nifti
		anat[ "anat" ] = os.path.join( anat[ "dir" ],
		                               "%s_anat" % subj_conf[ "subj_id" ]
		                             )

		return anat


	def get_roi_paths( subj_dir, roi_names ):
		""" Returns the paths for the ROIs

		Parameters
		----------
		subj_dir : string
			Path to subject base directory
		roi_names : list of strings
			Names of the ROIs

		Returns
		-------
		roi : dict of paths
			[ "dir" ] :
				Base directory for the ROI data.
			[ "mat" ] :
				Filenames for the ROI MAT files.
			[ "orig" ] :
				Filenames of the original ROI images.
			[ "rs" ] :
				Filenames of the resampled ROI images.

		"""

		roi = {}

		# base ROI directory
		roi[ "dir" ] = os.path.join( subj_dir, "roi" )

		# path to mat file for each ROI
		roi[ "mat" ] = [ os.path.join( roi[ "dir" ],
		                               "%s.mat" % roi_name
		                             )
		                 for roi_name in roi_names
		               ]

		# path to the nifti file of each ROI, in its original (anatomical) image
		# space
		roi[ "orig" ] = [ os.path.join( roi[ "dir" ],
		                               "%s_mask" % roi_name
		                              )
		                  for roi_name in roi_names
		                ]

		# path to the nifti file of each ROI in the same space as the functional
		# images
		roi[ "rs" ] = [ os.path.join( roi[ "dir" ],
		                              "rs_%s_mask" % roi_name
		                            )
		                for roi_name in roi_names
		              ]

		return roi


	analysis_conf = get_conf()[ "ana" ]

	subj_conf = get_subj_conf()[ subj_id ]

	paths = {}

	paths[ "exp" ] = get_exp_paths()

	paths[ "subj" ] = { "dir" : os.path.join( paths[ "exp" ][ "dir" ],
	                                          "subj_data",
	                                          subj_conf[ "subj_id" ]
	                                        )
	                  }

	paths[ "func" ] = get_func_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "fmap" ] = get_fmap_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "anat" ] = get_anat_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "roi" ] = get_roi_paths( paths[ "subj" ][ "dir" ],
	                                analysis_conf[ "rois" ]
	                              )

	paths[ "design" ] = get_design_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "ana" ] = get_analysis_paths( paths[ "subj" ][ "dir" ] )

	paths[ "task" ] = get_task_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	return paths

