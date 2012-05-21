"""Glass coherence block design fMRI experiment.
"""

from __future__ import division

import os

import numpy as np

import psychopy.visual, psychopy.filters, psychopy.misc, psychopy.event
import psychopy.core, psychopy.log

import psychopy.logging
psychopy.logging.console.setLevel( psychopy.logging.CRITICAL )

import glass_coherence_block.config
import stimuli.psychopy_ext, stimuli.utils

__author__ = "Damien Mannion"
__license__ = "GPL"
__version__ = "1"
__maintainer__ = "Damien Mannion"
__email__ = "dmannion@umn.edu"
__status__ = "Collection"


def run( subj_id, run_num ):
	"""Execute a run of the Glass coherence block design fMRI experiment.

	Parameters
	----------
	subj_id : string
		Subject ID, in the Olman lab system.
	run_num : integer,  1 <= run_num <= 12
		Run number.

	Returns
	-------
	status : integer, { 0, 1 }
		Returns 0 if run completed, 1 if aborted.

	"""

	# get the paths that will hold the run sequence and task information,
	# checking that they don't exist already
	( seq_path, task_path, resp_path ) = get_log_paths( subj_id, run_num )

	# load the experiment configuration
	conf = glass_coherence_block.config.get_conf()

	# get the event sequence info
	seq = get_seq( conf, run_num )

	# get the dictionary that tells us what the columns in the sequence mean
	seq_ind = get_seq_ind()

	# initialise the behavioural task details
	( task, targets ) = init_task( conf )
	task_resp = []

	# initialise the display window
	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = conf[ "acq" ][ "monitor_name" ],
	                              fullscr = True,
	                              allowGUI = False
	                            )

	try:
		stim = init_stim( conf, win )
	except:
		win.close()
		raise

	fix_text = psychopy.visual.TextStim( win = win,
	                                     text = "",
	                                     height = 26,
	                                     units = "pix",
	                                     bold = False
	                                   )

	target_pos = ( -40, 40 )

	target_text = [ psychopy.visual.TextStim( win = win,
	                                          text = "%s" % targets[ i, 0 ],
	                                          height = 26,
	                                          units = "pix",
	                                          pos = ( target_pos[ i ], 0 ),
	                                          color = np.repeat( targets[ i, 1 ], 3 )
	                                        )
	                for i in xrange( 2 )
	              ]

	_ = [ t_text.draw() for t_text in target_text ]


	# initialise the clock that keeps track of how long the run has been going
	run_clock = psychopy.core.Clock()

	# set the keys
	quit_key = 'q'
	trigger_key = 't'

	# wait for the trigger

	win.flip()

	k = psychopy.event.waitKeys( keyList = [ quit_key, trigger_key ] )

	# start the clock
	run_clock.reset()

	# quit, if the user wanted to
	if quit_key in k:
		print "User aborted"
		win.close()
		return 1

	run_time = run_clock.getTime()

	# keep looping until the time has elapsed
	while run_time < conf[ "exp" ][ "run_full_len_s" ]:

		i_evt = np.where( run_time > seq[ :, seq_ind[ "time_s" ] ] )[ 0 ]

		if i_evt.shape[ 0 ] > 0:

			i_evt = i_evt[ -1 ]

			evt_time = run_time - seq[ i_evt, seq_ind[ "time_s" ] ]

			if seq[ i_evt, seq_ind[ "generated" ] ] == 0:

				ori = seq[ i_evt, seq_ind[ "ori" ] ]
				coh = seq[ i_evt, seq_ind[ "coh" ] ]
				contrast = seq[ i_evt, seq_ind[ "contrast" ] ]

				stim.set_ori( ori )
				stim.set_coh( coh )
				stim.set_contrast( contrast )

				stim.instantiate()

				seq[ i_evt, seq_ind[ "generated" ] ] = 1

			# if the event time is less than how long we want the stimulus up for, draw
			# the stimuli
			if evt_time < conf[ "exp" ][ "evt_stim_s" ]:
				stim.draw()

		i_task_evt = np.where( run_time > task[ :, 0 ] )[ 0 ][ -1 ]

		fix_text.setText( str( int( task[ i_task_evt, 1 ] ) ) )
		fix_text.setColor( np.repeat( task[ i_task_evt, 2 ], 3 ) )

		fix_text.draw()

		# draw to the screen
		win.flip()

		# get any responses
		keys = psychopy.event.getKeys( timeStamped = run_clock )

		for ( key, timestamp ) in keys:

			if key == quit_key:
				print "User abort"
				win.close()
				return 1

			else:
				task_resp.extend( [ ( key, timestamp ) ] )

		run_time = run_clock.getTime()

	# all done, time for cleanup
	# first, close the window
	win.close()

	# now save the sequence
	np.save( seq_path, seq )
	np.save( task_path, task )

	# convert the response list into an np array for saving
	task_resp_np = np.array( task_resp,
	                         dtype = [ ( "key", "S10" ), ( "time", float ) ]
	                       )
	# and save
	np.save( resp_path, task_resp_np )

	return 0


def init_stim( conf, win ):
	"""Initialises the stimuli for the experiment.

	Parameters
	----------
	conf : dict
		Experiment configuration, as output from ``glass_coherence_block.config.get_conf()``
	win : Window object
		PsychoPy window.

	Returns
	-------
	stim : list of PatchStim-s
		Outer list is for a given image, inner list is for the image regions.

	"""

	stim_class = stimuli.psychopy_ext.GlassPattern

	stim = stim_class( win = win,
	                   diam = conf[ "stim" ][ "diam_deg" ],
	                   n_dipoles = conf[ "stim" ][ "n_dipoles" ],
	                   pole_sep = conf[ "stim" ][ "pole_sep_deg" ],
	                   dot_size = conf[ "stim" ][ "dot_size_deg" ],
	                   ori_type = "polar",
	                   ori_deg = 0,
	                   coh = 0,
	                   units = "deg",
	                   mask_type = "annulus",
	                   mask_in = conf[ "stim" ][ "mask_in_deg" ],
	                   mask_out = conf[ "stim" ][ "mask_out_deg" ],
	                   mask_fringe = conf[ "stim" ][ "mask_fringe_deg" ]
	                 )

	return stim


def get_seq_ind():
	"""Defines a lookup table for the columns in the experiment sequence array.

	Returns
	-------
	seq_ind : dict, containing items (all ints):
		time_s : event onset time, in seconds
		block_num : run block number
		block_type : whether an 'A' or 'B' block
		coh : pattern coherence for this event
		ori : pattern orientation for this event
		contrast : pattern contrast for this event
		generated : whether the stimulus has been instantiated for this event

	"""

	seq_ind = { "time_s" : 0,
	            "block_num" : 1,
	            "block_type" : 2,
	            "coh" : 3,
	            "ori" : 4,
	            "contrast" : 5,
	            "generated" : 6
	          }

	return seq_ind


def get_seq( conf, run_num ):
	"""Get a sequence of events that consitute a run.

	Returns
	-------
	seq : ( evt, parameter ) numpy array
		Details for each event. The columns are given by ``get_seq_ind``
	run_num : int, one-based
		Run number.

	"""

	block_seq = get_block_seq( run_num )

	seq_ind = get_seq_ind()

	# init the empty sequence
	seq = np.empty( ( conf[ "exp" ][ "n_evt_per_run" ],
	                  len( seq_ind )
	                )
	              )
	seq.fill( np.NAN )

	# loop through each event
	for i_evt in xrange( conf[ "exp" ][ "n_evt_per_run" ] ):

		# onset time
		time_s = i_evt * conf[ "exp" ][ "evt_len_s" ]

		# add the pre-period duration
		time_s += conf[ "exp" ][ "pre_len_s" ]

		# block index
		i_block = np.floor( i_evt *
		                    conf[ "exp" ][ "evt_len_s" ] /
		                    conf[ "exp" ][ "block_len_s" ]
		                  )

		block_type = block_seq[ i_block ]

		if block_type > 0:

			coh = conf[ "stim" ][ "coh_levels" ][ block_type - 1 ]

			oris = conf[ "stim" ][ "ori_deg" ].copy()
			np.random.shuffle( oris )
			ori = oris[ 0 ]

			contrast = 1

		else:
			coh = 0
			ori = -1
			contrast = 0

		seq[ i_evt, seq_ind[ "time_s" ] ] = time_s
		seq[ i_evt, seq_ind[ "block_num" ] ] = i_block + 1
		seq[ i_evt, seq_ind[ "block_type" ] ] = block_type
		seq[ i_evt, seq_ind[ "coh" ] ] = coh
		seq[ i_evt, seq_ind[ "ori" ] ] = ori
		seq[ i_evt, seq_ind[ "contrast" ] ] = contrast
		seq[ i_evt, seq_ind[ "generated" ] ] = 0

	# some sanity checks
	assert( np.sum( np.isnan( seq ) ) == 0 )

	return seq


def get_block_seq( run_num ):
	"""Counterbalanced block sequence for four conditions + baseline design.
	"""

	mini_blocks = np.array( ( [ 1, 2, 3, 4, 0 ],
	                          [ 4, 3, 2, 1, 0 ],
	                          [ 3, 1, 4, 2, 0 ],
	                          [ 2, 4, 1, 3, 0 ]
	                        )
	                      )

	i_perm = np.mod( run_num - 1, 4 )

	blk_seq = mini_blocks[ mini_blocks[ i_perm, :-1 ] - 1, : ].flatten()

	# add an extra blank at the start
	blk_seq = np.concatenate( ( [ 0 ], blk_seq ) )

	for i_type in xrange( 4 ):
		assert( sum( blk_seq == ( i_type + 1 ) ) == 4 )

	return blk_seq


def get_log_paths( subj_id, run_num ):
	"""Return the paths to run log files.

	Parameters
	----------
	subj_id : string
		Subject ID
	run_num : int (>0)
		Run number

	Returns
	-------
	seq_path : string
		Path to save the run image sequence.
	task_path : string
		Path to save the task and response sequence.

	Notes
	-----
	* If either of the files exist already, then an exception is thrown

	"""

	seq_path = os.path.join( "../logs",
	                         ( "%s_glass_coherence_block_seq_%d.npy" %
	                           ( subj_id, run_num )
	                         )
	                       )

	if os.path.exists( seq_path ):
		raise IOError( "Path %s already exists" % seq_path )

	task_path = os.path.join( "../logs",
	                          ( "%s_glass_coherence_block_task_%d.npy" %
	                            ( subj_id, run_num )
	                          )
	                        )

	if os.path.exists( task_path ):
		raise IOError( "Path %s already exists" % task_path )

	resp_path = os.path.join( "../logs",
	                          ( "%s_glass_coherence_block_resp_%d.npy" %
	                            ( subj_id, run_num )
	                          )
	                        )

	if os.path.exists( resp_path ):
		raise IOError( "Path %s already exists" % resp_path )

	return ( seq_path, task_path, resp_path )


def init_task( conf ):
	"""Initialises the task timing.

	Returns
	-------
	task_lut : numpy array, shape of ( evt x info )
		Task lookup table, where dim two is ( time_s, digit, polarity, target )
	targets : numpy array, shape of ( target, info )
		Target information, stored as ( digit, polarity )

	"""

	n_task_per_run = int( conf[ "exp" ][ "run_full_len_s" ] *
	                      conf[ "task" ][ "rate_hz" ]
	                    )

	task_set = conf[ "task" ][ "set" ]
	np.random.shuffle( task_set )

	targets = np.array( [ [ task_set[ i ], conf[ "task" ][ "polarity" ][ i ] ]
	                      for i in xrange( 2 )
	                    ]
	                  )

	# second dim is (time, digit, polarity, target or not)
	task_lut = np.empty( ( n_task_per_run, 4 ) )

	for i_evt in xrange( n_task_per_run ):

		time_s = i_evt * ( 1.0 / conf[ "task" ][ "rate_hz" ] )

		curr_task_set = task_set.copy()
		curr_task_set = curr_task_set[ curr_task_set != task_lut[ i_evt - 1, 1 ] ]

		digit = curr_task_set[ np.random.randint( len( curr_task_set ) ) ]

		polarity = conf[ "task" ][ "polarity" ][ np.random.randint( 2 ) ]

		if np.any( np.logical_and( targets[ :, 0 ] == digit,
		                           targets[ :, 1 ] == polarity
		                         )
		         ):
			target = 1
		else:
			target = 0

		task_lut[ i_evt, : ] = [ time_s, digit, polarity, target ]

	return ( task_lut, targets )
