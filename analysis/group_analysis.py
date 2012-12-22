"""
Set of routines to analyse group fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import os.path
import logging

import numpy as np
import scipy.stats

import fmri_tools.utils

import glass_coherence_block.config
import glass_coherence_block.analysis.paths


def task_anova( conf, paths ):
	"""Analyse the task performance across subjects."""

	logger = logging.getLogger( __name__ )
	logger.info( "Running task statistics ..." )

# A,B fixed; C random;  AxB,BxC,C(A)
	anova_type = "5"

	# stimulus conditions
	n_a_levels = len( conf[ "stim" ][ "coh_levels" ] )

	# time bin
	n_b_levels = conf[ "ana" ][ "task_perf_n_bins" ]

	# subjects
	n_c_levels = len( conf[ "all_subj" ] )

	# design spec
	cmd = [ "3dANOVA3",
	        "-DAFNI_FLOATIZE=YES",  # why not
	        "-overwrite",
	        "-type", anova_type,
	        "-alevels", "{n:d}".format( n = n_a_levels ),
	        "-blevels", "{n:d}".format( n = n_b_levels ),
	        "-clevels", "{n:d}".format( n = n_c_levels ),
	      ]

	bucket = paths.task_anova.full( ".niml.dset" )

	# analysis / output options
	cmd.extend( [ "-fa", "cond",  # main effect of condition
	              "-fb", "time",  # main effect of time bin
	              "-fab", "cond_x_time",  # interaction between condition and time bin
	              "-bucket", bucket
	            ]
	          )

	# data input
	subj_ids = conf[ "all_subj" ].keys()

	for i_c in xrange( n_c_levels ):

		subj_id = subj_ids[ i_c ]

		subj_conf = glass_coherence_block.config.get_conf( subj_id )
		subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

		for i_a in xrange( n_a_levels ):
			for i_b in xrange( n_b_levels ):

				dset_path = subj_paths.task.data.full( "_{a:d}.niml.dset[{b:d}]".format( a = i_a,
				                                                                         b = i_b
				                                                                       )
				                                     )

				cmd.extend( [ "-dset",
				              "{a:d}".format( a = ( i_a + 1 ) ),
				              "{b:d}".format( b = ( i_b + 1 ) ),
				              "{c:d}".format( c = ( i_c + 1 ) ),
				              dset_path
				            ]
				          )

	os.chdir( paths.base.full() )

	fmri_tools.utils.run_cmd( " ".join( cmd ) )

	# now, to print a summary
	summ = "Summary:\n"

	dof_summ = fmri_tools.utils.get_dset_dof( bucket )

	dof_summ = [ dof_summ[ i ] for i in [ 1, 3, 5 ] ]

	results_cmd = [ "3dmaskdump",
	                "-noijk",
	                "-nozero",
	                bucket
	              ]

	results = fmri_tools.utils.run_cmd( " ".join( results_cmd ) )
	results = results.std_out.strip().split( " " )

	results = [ results[ i ] for i in [ 1, 3, 5 ] ]

	dof = [ dof_summ[ i ].partition( "(" )[ -1 ].strip( ")" ).split( "," )
	        for i in [ 0, 1, 2 ]
	      ]

	p = [ fmri_tools.utils.get_f_p( float( results[ i ] ),
	                                float( dof[ i ][ 0 ] ),
	                                float( dof[ i ][ 1 ] )
	                              )
	      for i in xrange( 3 )
	    ]

	# main effect of condition
	summ += ( "\tMain effect of condition: " +
	          dof_summ[ 0 ] + " = " + results[ 0 ] + ", p = " + p[ 0 ] + "\n"
	        )

	# main effect of bin
	summ += ( "\tMain effect of bin: " +
	          dof_summ[ 1 ] + " = " + results[ 1 ] + ", p = " + p[ 1 ] + "\n"
	        )

	# condition x bin
	summ += ( "\tInteraction between condition and bin: " +
	          dof_summ[ 2 ] + " = " + results[ 2 ] + ", p = " + p[ 2 ] + "\n"
	        )

	fmri_tools.utils.write_to_log( summ )


def roi_prep( conf, paths ):
	"""Prepare the ROI values"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running ROI preparation..." )

	n_subj = len( conf[ "all_subj" ] )
	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	for ( roi_name, roi_id ) in conf[ "ana" ][ "rois" ]:

		logger.info( "\tConsidering {r:s}...".format( r = roi_name ) )

		psc_vals = np.empty( ( n_subj, n_cond ) )
		psc_vals.fill( np.NAN )

		for ( i_subj, subj_id ) in enumerate( conf[ "all_subj" ] ):

			subj_conf = glass_coherence_block.config.get_conf( subj_id )
			subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

			# nodes x ( roi id, .. )
			subj_psc = np.loadtxt( subj_paths.roi.psc.full( ".txt" ) )

			i_roi_nodes = ( subj_psc[ :, 0 ] == int( roi_id ) )

			subj_psc_mean = np.mean( subj_psc[ i_roi_nodes, 1: ], axis = 0 )

			psc_vals[ i_subj, : ] = subj_psc_mean

		assert( np.sum( np.isnan( psc_vals ) ) == 0 )

		psc_path = paths.psc.full( "_{roi:s}.txt".format( roi = roi_name ) )

		# save
		np.savetxt( psc_path, psc_vals )

		# now to normalise by subtracting the subject mean
		norm_psc_vals = np.empty( psc_vals.shape )
		norm_psc_vals.fill( np.NAN )

		# could use expansion, but just to be safe do it the verbose way
		for i_subj in xrange( n_subj ):

			subj_mean = np.mean( psc_vals[ i_subj, : ] )

			norm_psc_vals[ i_subj, : ] = ( psc_vals[ i_subj, : ] - subj_mean )

		assert( np.sum( np.isnan( norm_psc_vals ) ) == 0 )

		norm_psc_path = paths.psc.full( "_{roi:s}-norm.txt".format( roi = roi_name ) )

		np.savetxt( norm_psc_path, norm_psc_vals )



def roi_perm( conf, paths ):
	"""Permutation test on the ROI values"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running ROI permutations..." )

	n_subj = len( conf[ "all_subj" ] )
	n_cond = len( conf[ "stim" ][ "coh_levels" ] )

	n_con = len( conf[ "ana" ][ "con_coefs" ] )
	n_perm = conf[ "ana" ][ "n_perm" ]

	roi_seeds = conf[ "ana" ][ "roi_seeds" ]

	for ( ( roi_name, _ ), roi_seed ) in zip( conf[ "ana" ][ "rois" ], roi_seeds ):

		logger.info( "\tConsidering {r:s}...".format( r = roi_name ) )

		# might as well used normed; shouldn't matter
		psc = np.loadtxt( paths.psc.full( "_{roi:s}-norm.txt".format( roi = roi_name ) ) )

		con_data = np.empty( ( n_perm + 1, n_con ) )
		con_data.fill( np.NAN )

		np.random.seed( roi_seed )

		# first, handle the unpermuted data
		con_data[ 0, : ] = [ np.sum( np.mean( psc, axis = 0 ) * con_coef )
		                     for con_coef in conf[ "ana" ][ "con_coefs" ]
		                   ]

		# now the permuted
		for i_perm in xrange( n_perm ):

			perm_psc = [ psc[ i_subj, np.random.permutation( n_cond ) ]
			             for i_subj in xrange( n_subj )
			           ]

			con_data[ i_perm + 1, : ] = [ np.sum( np.mean( perm_psc, axis = 0 ) * con_coef )
			                              for con_coef in conf[ "ana" ][ "con_coefs" ]
			                            ]

		# now to save
		con_data_path = paths.con_data.full( "_{roi:s}.txt".format( roi = roi_name ) )

		np.savetxt( con_data_path, con_data )


def roi_stat( conf, paths ):
	"""Descriptive and inferential stats on the ROI data"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running ROI permutations..." )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		logger.info( "\tConsidering {r:s}...".format( r = roi_name ) )

		# descriptive stats
		psc = np.loadtxt( paths.psc.full( "_{roi:s}-norm.txt".format( roi = roi_name ) ) )

		descrip = np.empty( ( 2, psc.shape[ 1 ] ) )
		descrip.fill( np.NAN )

		# mean
		descrip[ 0, : ] = np.mean( psc, axis = 0 )
		# sem
		descrip[ 1, : ] = scipy.stats.sem( psc, axis = 0 )

		assert( np.sum( np.isnan( descrip ) ) == 0 )

		descrip_path = paths.descrip.full( "_{roi:s}.txt".format( roi = roi_name ) )

		np.savetxt( descrip_path, descrip )

		# inferential stats

		# ( 1 + n_perms ) x ( lin, quad, cub )
		con_data = np.loadtxt( paths.con_data.full( "_{roi:s}.txt".format( roi = roi_name ) ) )

		# this gives `p` as between 0 and 100
		p = np.array( [ scipy.stats.percentileofscore( a = con_data[ 1:, i_con ],
		                                               score = con_data[ 0, i_con ]
		                                             )
		                for i_con in xrange( con_data.shape[ 1 ] )
		              ]
		            )

		# flip around, convert to [ 0, 1 ], and double (for two-tailed)
		p = ( 100.0 - p ) / 100.0 * 2.0

		con_stat_path = paths.con_stat.full( "_{roi:s}.txt".format( roi = roi_name ) )

		np.savetxt( con_stat_path, p, "%.4f" )

