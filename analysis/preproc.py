"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path
import tempfile

import numpy as np

import fmri_tools.preproc, fmri_tools.utils


def convert( paths, conf ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# aggregate the dicom directories
	raw_dirs = ( paths[ "func" ][ "raw_dirs" ] +
	             paths[ "fmap" ][ "raw_mag_dirs" ] +
	             paths[ "fmap" ][ "raw_ph_dirs" ]
	           )

	# aggregate the output directories
	nii_dirs = ( paths[ "func" ][ "run_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ]
	           )

	# aggregate the images paths
	img_paths = ( paths[ "func" ][ "orig_files" ] +
	              paths[ "fmap" ][ "mag_files" ] +
	              paths[ "fmap" ][ "ph_files" ]
	            )

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


def st_motion_correct( paths, conf ):
	"""Performs slice-timing and motion correction"""

	# get the order of runs to pass to the correction algorithm
	# this is done because the algorithm realigns all to the first entry, which
	# normally corresponds to the first run acquired, but we might want them to
	# be aligned with a different run - one closer to a fieldmap, for example
	run_order = conf[ "subj" ][ "run_st_mot_order" ]

	# reorder the paths
	# (the -1 is because the runs are specified in subj_conf in a one-based
	# index; ie. run 1 is the first run)
	orig_paths = [ paths[ "func" ][ "orig_files" ][ i_run - 1 ]
	               for i_run in run_order
	             ]
	corr_paths = [ paths[ "func" ][ "corr_files" ][ i_run - 1 ]
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
	                                                   slice_info
	                                                 )

	# save the estimated motion parameters
	np.save( paths[ "summ" ][ "mot_est_file" ],
	         arr = motion_est
	       )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( corr_paths,
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

	# combine the experiment and localiser functional info
	func_fmap = paths[ "func" ][ "fmap_files" ]

	func_corr = paths[ "func" ][ "corr_files" ]

	func_uw = paths[ "func" ][ "uw_files" ]

	for i_run in xrange( len( func_corr )  ):

		fmri_tools.preproc.unwarp( func_corr[ i_run ],
		                           func_fmap[ i_run ],
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


def trim( paths, conf ):
	"""Trims the timecourses"""

	exp_start_vol = conf[ "exp" ][ "pre_len_s" ] / conf[ "acq" ][ "tr_s" ]
	exp_n_vol = conf[ "exp" ][ "run_len_s" ] / conf[ "acq" ][ "tr_s" ]

	for ( uw_file, trim_file ) in zip( paths[ "func" ][ "uw_files" ],
	                                   paths[ "func" ][ "trim_files" ]
	                                 ):

		exp_trim_cmd = [ "fslroi",
		                 uw_file,
		                 trim_file,
		                 "%d" % exp_start_vol,
		                 "%d" % exp_n_vol
		               ]

		fmri_tools.utils.run_cmd( exp_trim_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )


def surf_reg( paths, conf ):
	"""Coregisters an anatomical with the SUMA reference"""

	fmri_tools.preproc.surf_reg( paths[ "reg" ][ "rs_exp_anat" ],
	                             paths[ "reg" ][ "surf_anat" ],
	                             paths[ "summ" ][ "log_file" ]
	                           )


def vol_to_surf( paths, conf ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	start_dir = os.getcwd()

	vol_files = paths[ "func" ][ "trim_files" ]

	surf_files = paths[ "func" ][ "surf_files" ]

	for ( vol_file, surf_file ) in zip( vol_files, surf_files ):

		file_dir = os.path.split( vol_file )[ 0 ]

		os.chdir( file_dir )

		for hemi in [ "lh", "rh" ]:

			out = "%s_%s.niml.dset" % ( surf_file,
			                            hemi
			                          )

			surf_cmd = [ "3dVol2Surf",
			             "-spec", "%s%s.spec" % ( paths[ "reg" ][ "spec" ], hemi ),
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "15",
			             "-f_index", "nodes",
			             "-f_p1_fr", "0.0",
			             "-f_pn_fr", "0.0",
			             "-sv", paths[ "reg" ][ "reg" ],
			             "-grid_parent", "%s.nii" % vol_file,
			             "-out_niml", out,
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( surf_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	os.chdir( start_dir )


def design_prep( paths, conf ):
	"""Prepares the designs for GLM analysis"""

	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	# exp
	run_files = [ open( "%s%d.txt" % ( paths[ "ana" ][ "exp_time_files" ],
	                                   cond_num
	                                 ),
	                    "w"
	                  )
	              for cond_num in np.arange( 1, n_cond + 1 )
	            ]

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		run_times = []
		run_conds = []

		run_seq = np.load( "%s%d.npy" % ( paths[ "log" ][ "seq_base" ], i_run + 1 ) )

		for i_evt in xrange( run_seq.shape[ 0 ] ):

			is_transition = ( run_seq[ i_evt, 1 ] != run_seq[ i_evt - 1, 1 ] )

			if is_transition:

				start_time_s = run_seq[ i_evt, 0 ]

				cond = int( run_seq[ i_evt, 2 ] )

				run_times.append( start_time_s )
				run_conds.append( cond )

		run_times = np.array( run_times )
		run_conds = np.array( run_conds )

		run_times -= conf[ "exp" ][ "pre_len_s" ]

		ok = np.logical_and( run_times >= 0,
		                     run_times < conf[ "exp" ][ "run_len_s" ]
		                   )

		# exclude fixation blocks
		ok = np.logical_and( ok, run_conds > 0 )

		run_times = run_times[ ok ]
		run_conds = run_conds[ ok ]

		for ( i_evt, run_time ) in enumerate( run_times ):

			run_files[ run_conds[ i_evt ] - 1 ].write( "%.5f\t" % run_time )

		_ = [ run_file.write( "\n" ) for run_file in run_files ]

	_ = [ run_file.close() for run_file in run_files ]
