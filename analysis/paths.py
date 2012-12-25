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
	paths.roi = _get_roi_paths( conf, paths )
	paths.task = _get_task_paths( conf, paths )

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

	ana.bltc = ana.base + ( file_base + "bltc" )
	ana.bl = ana.base + ( file_base + "bl" )
	ana.psc = ana.base + ( file_base + "psc" )

	ana.fdr = ana.base + ( file_base + "fdr" )
	ana.mask = ana.base + ( file_base + "mask" )

	return ana


def _get_roi_paths( conf, paths ):
	"""Get the paths for the ROI analysis"""

	roi = fmri_tools.paths.PathsHandler()

	roi.base = paths.base / "rois"

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id, exp_id = exp_id )

	roi.vl = roi.base + ( file_base + "vis_loc_rois" )
	roi.vl_subset = roi.base + ( file_base + "vis_loc_rois_subset" )

	roi.mask_rois = roi.base + ( file_base + "mask_rois" )
	roi.rois = roi.base + ( file_base + "rois" )
	roi.dv_rois = roi.base + ( file_base + "dorsal_ventral_rois" )

	wedge_file = "{s:s}_vis_loc_{d:s}_wedge-angle".format( s = subj_id,
	                                                       d = conf[ "subj" ][ "wedge_date" ]
	                                                     )

	roi.fs_wedge = fmri_tools.paths.Path( os.path.join( "/labs/olmanlab/FsVisLoc/",
	                                                    subj_id,
	                                                    conf[ "subj" ][ "wedge_date" ],
	                                                    "dt", "wedge"
	                                                  ),
	                                      wedge_file
	                                    )

	roi.wedge = roi.base + ( file_base + "wedge" )

	roi.psc = roi.base + ( file_base + "psc" )

	roi.vf = roi.base + ( file_base + "vf" )

	roi.con_coef = roi.base + ( file_base + "con_coef" )

	roi.lin_v3_vf = roi.base + ( file_base + "lin_v3_vf" )

	return roi


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


def _get_task_paths( conf, paths ):
	"""Get the paths for the task analysis"""

	task = fmri_tools.paths.PathsHandler()

	task.base = paths.base / "task"

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	file_base = "{subj_id:s}_{exp_id:s}-task_".format( subj_id = subj_id, exp_id = exp_id )

	task.perf = task.base + ( file_base + "perf" )
	task.data = task.base + ( file_base + "data" )

	return task


def get_group_paths( conf ):
	"""Get the path structure for the group analysis"""

	grp = fmri_tools.paths.PathsHandler()

	grp.base = fmri_tools.paths.Path( "/labs/olmanlab/Data7T/GlassCoherenceBlock/group_data" )

	file_base = "glass_coherence_block"

	grp.log = grp.base + ( file_base + "-log.log" )

	grp.task_anova = grp.base + ( file_base + "-task" )

	grp.psc = grp.base + ( file_base + "-psc" )

	grp.con_data = grp.base + ( file_base + "-con_data" )

	grp.descrip = grp.base + ( file_base + "-descrip" )

	grp.con_stat = grp.base + ( file_base + "-con_stat" )

	# figures

	grp.fig_base = grp.base / "figures"

	grp.fig_task = grp.fig_base + ( file_base + "-task" )
	grp.fig_psc = grp.fig_base + ( file_base + "-psc" )
	grp.fig_roi_psc = grp.fig_base + ( file_base + "-roi_psc" )

	return grp
