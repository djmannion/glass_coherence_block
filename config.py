"""Configuration for the Glass coherence block design fMRI experiment.
"""

from __future__ import division

import numpy as np

import fmri_tools.utils


def get_conf( subj_id = None ):
	"""Overall experiment configuration.

	Returns
	-------
	conf : dict, with items:
		exp : holds overall experiment configurations
		stim : stimulus configuration
		task : behavioural task configuration
		acq : acquisition configuration

	"""

	conf = {}

	conf[ "exp" ] = _get_exp_conf()
	conf[ "stim" ] = _get_stim_conf()
	conf[ "task" ] = _get_task_conf()
	conf[ "acq" ] = _get_acq_conf()
	conf[ "ana" ] = _get_ana_conf( conf )
	conf[ "all_subj" ] = _get_subj_conf()

	if subj_id is not None:

		conf[ "subj" ] = _get_subj_conf( subj_id )

	return conf


def _get_ana_conf( conf ):
	"""
	"""

	rois = [ [ "v1", "1" ],
	         [ "v2", "2" ],
	         [ "v3", "3" ],
	         [ "hv4", "4" ],
	         [ "vo1", "10" ],
	         [ "v3ab", "5" ],
	         [ "lo1", "6" ],
	         [ "lo2", "7" ],
	         [ "hmtp", "9" ]
	      ]

	hrf_model = "SPMG1(%d)" % conf[ "exp" ][ "block_len_s" ]

	poly_ord = "3"

	loc_q = 0.01

	ana_conf = { "rois" : rois,
	             "poly_ord" : poly_ord,
	             "loc_q" : loc_q,
	             "hrf_model" : hrf_model
	           }

	return ana_conf


def _get_stim_conf():
	"""Get the stimulus configuration."""

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
	"""Gets the experiment configuration."""

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

	exp_conf[ "id" ] = "glass_coherence_block"

	return exp_conf


def _get_task_conf():
	"""Gets the task configuration."""

	task_conf = {}

	task_conf[ "set" ] = np.arange( 10 )

	task_conf[ "polarity" ] = ( -1, +1 )

	task_conf[ "rate_hz" ] = 3.0

	return task_conf


def _get_acq_conf():
	"""Get the acquisition configuration."""

	acq_conf = {}

	acq_conf[ "monitor_name" ] = "UMN_7T"

	acq_conf[ "tr_s" ] = 2.0

	acq_conf[ "delta_te_ms" ] = 1.02

	acq_conf[ "ras" ] = ( "-x", "-z", "-y" )

	acq_conf[ "dwell_ms" ] = 0.325

	acq_conf[ "slice_order" ] = fmri_tools.utils.get_slice_order( 36,
	                                                              "interleaved"
	                                                            )

	acq_conf[ "slice_axis" ] = 1

	acq_conf[ "slice_acq_dir" ] = "+1"

	acq_conf[ "ph_encode_dir" ] = "z"


	return acq_conf


def _get_subj_conf( subj_id = None ):
	"""Gets the configuration info for each subject."""

	s1021 = { "subj_id" : "s1021",
	          "acq_date" : "20120522",
	          "comments" : "",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
	          "node_k" : { "lh" : 140847,
	                       "rh" : 141381
	                     }
	        }

	s1011 = { "subj_id" : "s1011",
	          "acq_date" : "20120601",
	          "comments" : "",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
	          "node_k" : { "lh" : 128434,
	                       "rh" : 128461
	                     }
	        }

	s1000 = { "subj_id" : "s1000",
	          "acq_date" : "20120605",
	          "comments" : "",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
	          "node_k" : { "lh" : 130318,
	                       "rh" : 131151
	                     }
	        }

	s1008 = { "subj_id" : "s1008",
	          "acq_date" : "20120605",
	          "comments" : "",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
	          "node_k" : { "lh" : 140427,
	                       "rh" : 141898
	                     }
	        }

	s1010 = { "subj_id" : "s1010",
	          "acq_date" : "20120605",
	          "comments" : """Re-ran third run after subject reported falling
	                       asleep""",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
	          "node_k" : { "lh" : 124787,
	                       "rh" : 127339
	                     }
	        }

	s1032 = { "subj_id" : "s1032",
	          "acq_date" : "20120612",
	          "comments" : "",
	          "n_runs" : 12,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
	          "node_k" : { "lh" : 126078,
	                       "rh" : 126201
	                     }
	        }


	subj_conf = { "s1021" : s1021,
	              "s1000" : s1000,
	              "s1008" : s1008,
	              "s1010" : s1010,
	              "s1032" : s1032,
#	              "s1011" : s1011
	            }

	if subj_id is None:
		return subj_conf
	else:
		return subj_conf[ subj_id ]
