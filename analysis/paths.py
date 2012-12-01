"""Paths for the Glass coherence block design fMRI experiment"""

from __future__ import division

import os

import fmri_tools.paths


def get_subj_paths( conf ):
	"""Get the path structure for a given subject"""

	base_dir = os.path.join( "/labs/olmanlab/Data7T/GlassCoherenceBlock/subj_data/",
	                         conf[ "subj" ][ "subj_id" ]
	                       )

	paths = fmri_tools.paths.get_default_paths( base_dir = base_dir,
	                                            subj_id = conf[ "subj" ][ "subj_id" ],
	                                            study_id = conf[ "exp" ][ "id" ],
	                                            n_runs = conf[ "subj" ][ "n_runs" ]
	                                          )

	paths.ana = _get_ana_paths( conf, paths )
	paths.logs = _get_log_paths( conf, paths )

	return paths


def _get_ana_paths( conf, paths ):
	"""Get the paths for the analysis"""

	ana = fmri_tools.paths.PathsHandler()

	ana.base = paths.base / "analysis"

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id, exp_id = exp_id )

	ana.stim_times = ana.base + ( file_base + "stim_times" )

	ana.glm = ana.base + ( file_base + "glm" )
	ana.beta = ana.base + ( file_base + "beta" )

	return ana

def _get_log_paths( conf, paths ):
	"""Get the paths for the logfiles"""

	logs = fmri_tools.paths.PathsHandler()

	logs.base = paths.base / "logs"

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	file_base = "{subj_id:s}_{exp_id:s}_".format( subj_id = subj_id, exp_id = exp_id )

	logs.seq = logs.base + ( file_base + "seq" )
	logs.task = logs.base + ( file_base + "task" )
	logs.resp = logs.base + ( file_base + "resp" )

	return logs


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
