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
