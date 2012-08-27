"""Paths for the Glass coherence block design fMRI experiment"""

from __future__ import division

import os

import fmri_tools.paths


def _get_study_paths():
	"""Get the path structure for the study"""

	base_dir = "/labs/olmanlab/Data7T/GlassCoherenceBlock/"

	study_paths = fmri_tools.paths.get_study_paths( base_dir )

	return study_paths


def _get_func_paths( conf, paths ):
	"Returns the paths associated with the functional acquisitions"

	func_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         conf[ "subj" ][ "subj_id" ],
	                         "func"
	                       )

	func_exp_dir = os.path.join( func_dir )

	func_paths = fmri_tools.paths.get_func_paths( func_exp_dir,
	                                              conf[ "subj" ][ "subj_id" ],
	                                              conf[ "subj" ][ "n_runs" ],
	                                              conf[ "exp" ][ "id" ]
	                                            )

	func_paths[ "mask_files" ] = [ orig_file.replace( "orig", "mask" )
	                               for orig_file in func_paths[ "orig_files" ]
	                             ]

	func_paths[ "ivar_files" ] = [ orig_file.replace( "orig", "ivar" )
	                               for orig_file in func_paths[ "orig_files" ]
	                             ]

	func_paths[ "trim_files" ] = [ orig_file.replace( "orig", "trim" )
	                               for orig_file in func_paths[ "orig_files" ]
	                             ]

	func_paths[ "surf_files" ] = [ orig_file.replace( "orig", "surf" )
	                               for orig_file in func_paths[ "orig_files" ]
	                             ]

	paths[ "func" ] = func_paths

	return paths


def _get_summ_paths( conf, paths ):
	"Returns the paths for images summarising the functional session"

	summ_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         conf[ "subj" ][ "subj_id" ],
	                         "func"
	                       )

	summ_paths = fmri_tools.paths.get_func_summ_paths( summ_dir,
	                                                   conf[ "subj" ][ "subj_id" ],
	                                                   conf[ "exp" ][ "id" ]
	                                                 )

	summ_paths[ "mot_est_file" ] = summ_paths[ "mot_est_file" ].replace( "npy",
	                                                                     "txt"
	                                                                   )

	log_file = "%s_%s-log.log" % ( conf[ "subj" ][ "subj_id" ],
	                               conf[ "exp" ][ "id" ]
	                             )

	summ_paths[ "log_file" ] = os.path.join( paths[ "study" ][ "subj_dir" ],
	                                         conf[ "subj" ][ "subj_id" ],
	                                         log_file
	                                       )

	paths[ "summ" ] = summ_paths

	return paths


def _get_fmap_paths( conf, paths ):
	"Returns the paths for the session's fieldmaps"

	fmap_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         conf[ "subj" ][ "subj_id" ],
	                         "fmap"
	                       )

	fmap_paths = fmri_tools.paths.get_fmap_paths( fmap_dir,
	                                              conf[ "subj" ][ "subj_id" ],
	                                              conf[ "exp" ][ "id" ],
	                                              conf[ "subj" ][ "n_fmaps" ]
	                                            )

	paths[ "fmap" ] = fmap_paths

	return paths


def _get_reg_paths( conf, paths ):
	"Returns the paths for the session's registration with FreeSurfer anatomy"

	reg_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        conf[ "subj" ][ "subj_id" ],
	                        "reg"
	                      )

	reg_paths = {}

	subj_id = conf[ "subj" ][ "subj_id" ]

	reg_paths[ "anat" ] = os.path.join( reg_dir,
	                                    "%s_anat+orig" % subj_id
	                                  )

	reg_paths[ "reg_anat" ] = os.path.join( reg_dir,
	                                        "%s_reg_anat+orig" % subj_id
	                                      )

	reg_paths[ "rs_exp_anat" ] = os.path.join( reg_dir,
	                                           "%s_anat_rs+orig" % subj_id
	                                         )

	surf_dir = os.path.join( os.environ[ "SUBJECTS_DIR" ],
	                         subj_id,
	                         "SUMA"
	                       )

	reg_paths[ "surf_anat" ] = os.path.join( surf_dir,
	                                         "%s_SurfVol+orig" % subj_id
	                                       )

	reg_paths[ "spec" ] = os.path.join( surf_dir,
	                                    "%s_" % subj_id
	                                  )

	paths[ "reg" ] = reg_paths

	return paths


def _get_log_paths( conf, paths ):
	"""Get the paths for the logfiles"""

	subj_id = conf[ "subj" ][ "subj_id" ]

	log = {}

	log_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        subj_id,
	                        "logs"
	                      )

	log[ "base_dir" ] = log_dir

	log[ "seq_base" ] = os.path.join( log_dir,
	                                  "%s_glass_coherence_block_seq_" % subj_id
	                                )

	log[ "task_base" ] = os.path.join( log_dir,
	                                   "%s_glass_coherence_block_task_" % subj_id
	                                 )

	paths[ "log" ] = log

	return paths


def _get_ana_paths( conf, paths ):
	"""Get the paths for the analysis"""

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	ana = {}

	ana[ "base_dir" ] = os.path.join( paths[ "study" ][ "subj_dir" ],
	                                  subj_id,
	                                  "analysis"
	                                )

	ana[ "time_files" ] = os.path.join( ana[ "base_dir" ],
	                                    "%s_%s_stim_times-" % ( subj_id, exp_id )
	                                  )

	ana[ "mot_est" ] = os.path.join( ana[ "base_dir" ],
	                                 "%s_%s-mot_est.txt" % ( subj_id, exp_id )
	                               )

	ana[ "bl_poly" ] = os.path.join( ana[ "base_dir" ],
	                                 "%s_%s-bl_poly.txt" % ( subj_id, exp_id )
	                               )

	ana[ "pred_adj" ] = os.path.join( ana[ "base_dir" ],
	                                  "%s_%s-pred_adj" % ( subj_id, exp_id )
	                                )

	ana[ "glm" ] = os.path.join( ana[ "base_dir" ],
	                             "%s_%s-glm" % ( subj_id, exp_id )
	                           )

	ana[ "beta" ] = os.path.join( ana[ "base_dir" ],
	                              "%s_%s-beta" % ( subj_id, exp_id )
	                            )

	ana[ "bltc" ] = os.path.join( ana[ "base_dir" ],
	                              "%s_%s-bltc" % ( subj_id, exp_id )
	                            )

	ana[ "bl" ] = os.path.join( ana[ "base_dir" ],
	                            "%s_%s-bl" % ( subj_id, exp_id )
	                          )

	ana[ "psc" ] = os.path.join( ana[ "base_dir" ],
	                             "%s_%s-psc" % ( subj_id, exp_id )
	                           )

	ana[ "loc_fdr" ] = os.path.join( ana[ "base_dir" ],
	                                 "%s_%s-loc_fdr" % ( subj_id, exp_id )
	                               )

	ana[ "loc_mask" ] = os.path.join( ana[ "base_dir" ],
	                                  "%s_%s-loc_mask" % ( subj_id, exp_id )
	                                )

	ana[ "raw" ] = os.path.join( ana[ "base_dir" ],
	                             "%s_%s-raw" % ( subj_id, exp_id )
	                           )

	ana[ "raw_adj" ] = os.path.join( ana[ "base_dir" ],
	                                 "%s_%s-raw_adj" % ( subj_id, exp_id )
	                               )

	paths[ "ana" ] = ana

	return paths


def _get_roi_paths( conf, paths ):
	"""Get the paths for the ROI data"""

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	rois = {}

	rois[ "base_dir" ] = os.path.join( paths[ "study" ][ "subj_dir" ],
	                                   subj_id,
	                                   "rois"
	                                 )

	rois[ "dset" ] = os.path.join( rois[ "base_dir" ],
	                               "%s_%s-rois" % ( subj_id, exp_id )
	                             )

	rois[ "psc" ] = os.path.join( rois[ "base_dir" ],
	                              "%s_%s-psc" % ( subj_id, exp_id )
	                            )

	rois[ "raw_adj_tc" ] = os.path.join( rois[ "base_dir" ],
	                                     "%s_%s-raw_adj_tc" % ( subj_id, exp_id )
	                                   )

	rois[ "pred_adj_tc" ] = os.path.join( rois[ "base_dir" ],
	                                      "%s_%s-pred_adj_tc" % ( subj_id, exp_id )
	                                    )

	paths[ "rois" ] = rois

	return paths


def get_subj_paths( conf ):
	"""Get the path structure for a given subject"""

	paths = {}

	paths[ "study" ] = _get_study_paths()

	paths = _get_func_paths( conf, paths )

	paths = _get_summ_paths( conf, paths )

	paths = _get_fmap_paths( conf, paths )

	paths = _get_reg_paths( conf, paths )

	paths = _get_log_paths( conf, paths )

	paths = _get_ana_paths( conf, paths )

	paths = _get_roi_paths( conf, paths )

	return paths


def get_group_paths( conf ):
	"""Get the path structure for the group analysis"""

	paths = {}

	paths[ "study" ] = _get_study_paths()

	paths[ "base_dir" ] = os.path.join( paths[ "study" ][ "base_dir" ],
	                                    "group_data"
	                                  )

	paths[ "roi_mean" ] = os.path.join( paths[ "base_dir" ],
	                                    "glass_coherence_block_group"
	                                  )

	paths[ "roi_stat" ] = os.path.join( paths[ "base_dir" ],
	                                    "glass_coherence_block_group-stat"
	                                  )

	return paths
