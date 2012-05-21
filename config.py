"""Configuration for the Glass coherence block design fMRI experiment.
"""

from __future__ import division

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
	         "acq" : _get_acq_conf()
	       }

	return conf


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

#	s1000 = { "subj_id" : "s1000",
#	          "acq_date" : "20120106",
#	          "comments" : ""
#	        }

#	subj_conf = { "s1000" : s1000,
#	            }

	subj_conf = {}

	return subj_conf
