"""
Set of routines to analyse group fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import os.path
import tempfile

import numpy as np
import scipy.stats

import fmri_tools.utils

import glass_coherence_block.config
import glass_coherence_block.analysis.paths


def task_anova( conf, paths ):
	"""Analyse the task performance across subjects."""

# A,B fixed; C random;  AxB,BxC,C(A)
	anova_type = "5"

	# stimulus condition; coherences, plus 1 blank
	n_a_levels = len( conf[ "stim" ][ "coh_levels" ] ) + 1

	# time bin
	n_b_levels = conf[ "ana" ][ "task_perf_n_bins" ]

	# subjects
	n_c_levels = len( conf[ "all_subj" ] )

	# design spec
	cmd = [ "3dANOVA3",
	        "-DAFNI_FLOATIZE=YES",  # why not
	        "-overwrite",
	        "-type", "5",  # B fixed; C random;  AxB,BxC,C(A)
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


def roi_mean( paths, conf ):
	"""Average the nodes in each subjects ROIs"""

	roi_names = [ roi_info for ( roi_info, _ ) in conf[ "ana" ][ "rois" ] ]

	for roi_name in roi_names:

		# roi summary file, to write; one row per subject
		roi_path = "%s-%s.txt" % ( paths[ "roi_mean" ], roi_name )

		roi_file = open( roi_path, "w+" )

		for subj_id in conf[ "all_subj" ]:

			subj_conf = glass_coherence_block.config.get_conf( subj_id )
			subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

			roi_data = []

			for hemi in [ "lh", "rh" ]:

				# the file containing the psc values for this subject, ROI, and hemi
				roi_psc_file = "%s_%s_%s.txt" % ( subj_paths[ "rois" ][ "psc" ],
				                                  roi_name,
				                                  hemi
				                                )

				roi_data.append( np.loadtxt( roi_psc_file ) )

			# concatenate the two hemisphere lists over nodes
			roi_data = np.vstack( roi_data )

			# average over nodes
			roi_data = np.mean( roi_data, axis = 0 )

			# subtract the average over conditions
			roi_data -= np.mean( roi_data )

			# this is saving at the same precision as savetxt
			subj_roi_txt = ( subj_id +
			                 "\t%.18e\t%.18e\t%.18e\t%.18e\n" % tuple( roi_data )
			               )

			roi_file.write( subj_roi_txt )

		roi_file.close()


def roi_stat( paths, conf ):
	"""Permutation test on the ROI values"""

	roi_names = [ roi_info for ( roi_info, _ ) in conf[ "ana" ][ "rois" ] ]

	n_rois = len( roi_names )

	trends = np.array( [ [ -3, -1, +1, +3 ],  # linear
	                     [ +1, -1, -1, +1 ],  # quadratic
	                     [ -1, +3, -3, +1 ]   # cubic
	                   ]
	                 )

	( n_trends, _ ) = trends.shape

	# just check that they all sum to zero, like they should
	assert( np.all( np.sum( trends, axis = 1 ) == 0 ) )

	n_perm = conf[ "ana" ][ "n_perm" ]

	# roi x ( mean, sem, p, q ) x trend
	stat_data = np.empty( ( n_rois, 4, n_trends ) )

	roi_info = zip( roi_names, conf[ "ana" ][ "roi_seeds" ] )

	# iterate over each ROI
	for ( i_roi, ( roi_name, roi_seed ) ) in enumerate( roi_info ):

		roi_path = "%s-%s.txt" % ( paths[ "roi_mean" ], roi_name )

		# first column is subject id, so don't load it
		roi_data = np.loadtxt( roi_path, usecols = [ 1, 2, 3, 4 ] )

		( n_subj, n_cond ) = roi_data.shape

		# set the seed of the random number generator to a known value; subsequent
		# random calls will be reproducible
		np.random.seed( roi_seed )

		# initialise the permutation data container
		perm_trend = np.empty( ( n_perm, n_trends ) )
		perm_trend.fill( np.NAN )

		# iterate over each permutation
		for i_perm in xrange( n_perm ):

			# for each subject's data, randomly permute the condition indices
			# (without replacement)
			perm_data = [ roi_data[ i_subj, np.random.permutation( n_cond ) ]
			              for i_subj in xrange( n_subj )
			            ]

			# put in a useful form
			perm_data = np.vstack( perm_data )

			for ( i_trend, trend_coef ) in enumerate( trends ):

				# multiply by the trend coefficients, then sum over conditions
				trend_data = np.sum( perm_data * trend_coef, axis = 1 )

				# average over subjects as the statistic
				perm_trend[ i_perm, i_trend ] = np.mean( trend_data )

		# now we have our permutation distribution, we can calculate a p-value for
		# the measured values

		# first we need to know the measured trend values
		for ( i_trend, trend_coef ) in enumerate( trends ):

			# as above
			trend_data = np.sum( roi_data * trend_coef, axis = 1 )

			stat_data[ i_roi, 0, i_trend ] = np.mean( trend_data )

			stat_data[ i_roi, 1, i_trend ] = scipy.stats.sem( trend_data )

			# use the absolute value to get a two-tailed p
			p = scipy.stats.percentileofscore( np.abs( perm_trend[ :, i_trend ] ),
			                                   np.abs( stat_data[ i_roi, 0, i_trend ] )
			                                 )

			stat_data[ i_roi, 2, i_trend ] = ( 100 - p ) / 100.0

	# now we want to convert all the p's into q's, over ROIs
	for i_trend in xrange( n_trends ):

		p_file = tempfile.NamedTemporaryFile()
		q_file = tempfile.NamedTemporaryFile()

		# save the p data for all the ROIs for this trend
		np.savetxt( p_file.name, stat_data[ :, 2, i_trend ] )

		# use AFNI's '3dFDR' to compute q values
		fdr_cmd = [ "3dFDR",
		            "-input1D", p_file.name,
		            "-output", q_file.name,
		            "-qval",
		            "-new",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env()
		                        )

		# load the outputted q values
		q_vals = np.loadtxt( q_file.name )

		stat_data[ :, 3, i_trend ] = q_vals

	for ( i_roi, roi_name ) in enumerate( roi_names ):

		np.savetxt( "%s-%s.txt" % ( paths[ "roi_stat" ], roi_name ),
		            stat_data[ i_roi, :, : ]
		          )

