"""Plots figures for the Glass coherence block design fMRI experiment
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats

import glass_coherence_block.analysis.paths


def write_mask_cmap( rois, cmap_path ):
	"""Writes a SUMA colourmap for the mask display"""

	import brewer2mpl

	map_name = "Dark2"
	map_type = "qualitative"

	n_cols = len( rois )

	cols = brewer2mpl.get_map( "Dark2", "qualitative", n_cols ).mpl_colors

	i_rois = [ int( roi[ 1 ] ) for roi in rois ]

	cmap = np.zeros( ( max( i_rois ), 3 ) )

	for ( i, i_roi ) in enumerate( i_rois ):

		cmap[ i_roi - 1, : ] = cols[ i ]

	np.savetxt( cmap_path, cmap )


def plot_task_perf( conf, paths, save_path = None ):
	"""Plot the task performance"""

	cond_labels = [ "0", "33", "66", "100" ]

	n_cond = len( cond_labels )

	n_bins = conf[ "ana" ][ "task_perf_n_bins" ]

	time_bins = np.arange( start = 0,
	                       stop = conf[ "ana" ][ "task_bin_res_s" ] * n_bins,
	                       step = conf[ "ana" ][ "task_bin_res_s" ]
	                     )

	bin_offsets = np.array( [ -0.06, -0.02, 0.02, 0.06 ] ) / 2.5

	subj_ids = conf[ "all_subj" ].keys()

	n_subj = len( subj_ids )

	data = np.empty( ( n_subj, n_bins, n_cond ) )
	data.fill( np.NAN )

	for ( i_subj, subj_id ) in enumerate( subj_ids ):

		subj_conf = glass_coherence_block.config.get_conf( subj_id )
		subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

		data[ i_subj, ... ] = np.loadtxt( subj_paths.task.perf.full( ".txt" ) )

	# average over subjects; mean is 20 x 5
	mean = np.mean( data, axis = 0 )

	# standard error
	std_error = np.std( data, axis = 0 ) / np.sqrt( n_subj )

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 3.34646, 3.34646 * 0.75, forward = False )

	ax = fig.gca()
	ax.hold( True )

	for i_cond in xrange( n_cond ):

		x = time_bins + bin_offsets[ i_cond ]

		_ = [ ax.plot( [ x[ i ] ] * 2, [ mean[ i, i_cond ] - std_error[ i, i_cond ],
		                                 mean[ i, i_cond ] + std_error[ i, i_cond ]
		                               ],
		               c = [ 0.5 ] * 3
		             )
		      for i in xrange( mean.shape[ 0 ] )
		    ]

		ax.scatter( x,
		            mean[ :, i_cond ],
		            marker = "s",
		            edgecolor = [ 0 ] * 3,
		            facecolor = [ 0 ] * 4,
		            s = 5,
		            zorder = 100
		          )

	_cleanup_fig( ax )

	ax.set_xlim( ( -0.05, time_bins[ -1 ] + 0.05 ) )

	ax.set_xlabel( "Time from target onset (s)" )
	ax.set_ylabel( "Correlation (r)" )

	fig.tight_layout( pad = 0.5 )

	if save_path:
		fig.savefig( save_path )
	else:
		fig.show()



def _set_defaults():
	"""Set some sane defaults for figures.
	"""

	params = { 'axes.labelsize': 9 * ( 1 / 1.25 ),
	           'axes.titlesize' : 10,
	           'font.family' : 'Arial',
	           'font.sans-serif' : 'Helvetica',
	           'text.fontsize': 12,
	           'legend.fontsize': 7,
	           'xtick.labelsize': 8 * ( 1 / 1.25 ),
	           'xtick.direction' : 'out',
	           'xtick.major.size' : 2,
	           'ytick.labelsize': 8 * ( 1 / 1.25 ),
	           'ytick.direction' : 'out',
	           'ytick.major.size' : 2
	         }
	
	plt.rcParams.update( params )

	plt.ioff()


def _cleanup_fig( ax ):
	"""Apply some standard commands to clean up the axes on figures.
	"""

	for loc, spine in ax.spines.iteritems():

		spine.set_linewidth( 0.5 )

		if loc in [ "left", "bottom" ]:
			spine.set_position( ( "outward", 5 ) )
		elif loc in [ "right", "top" ]:
			spine.set_color( "none" )
		else:
			raise ValueError( "Unknown spine location: %s" % loc )

	ax.xaxis.set_ticks_position( "bottom" )
	ax.yaxis.set_ticks_position( "left" )


def plot_subj_tc( paths, conf, save_path = None ):
	"""Plots the raw and predicted (both adjusted) timecourses for each run and
	ROI for a single subject"""

	_set_defaults()

	n_vols_per_run = ( conf[ "exp" ][ "run_len_s" ] /
	                   conf[ "acq" ][ "tr_s" ]
	                 )

	y_head = 10

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		fig = plt.figure()

		fig.set_size_inches( 7, 7, forward = True )

		# assuming 12 runs
		gs = gridspec.GridSpec( 4, 3 )

		raw_data = [ np.loadtxt( "%s_%s_%s.txt" % ( paths[ "rois" ][ "raw_adj_tc" ],
		                                            roi_name,
		                                            hemi
		                                          )
		                       )
		             for hemi in [ "lh", "rh" ]
		           ]

		pred_data = [ np.loadtxt( "%s_%s_%s.txt" % ( paths[ "rois" ][ "pred_adj_tc" ],
		                                             roi_name,
		                                             hemi
		                                           )
		                        )
		             for hemi in [ "lh", "rh" ]
		           ]

		raw_data = np.vstack( raw_data )
		pred_data = np.vstack( pred_data )

		raw_data = np.mean( raw_data, axis = 0 )
		pred_data = np.mean( pred_data, axis = 0 )

		y_min = 0
		y_max = 0

		ax_list = []

		for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

			vol_range = np.arange( i_run * n_vols_per_run,
			                       i_run * n_vols_per_run + n_vols_per_run
			                     ).astype( "int" )

			ax = plt.subplot( gs[ i_run ] )

			ax.hold( True )

			ax.plot( raw_data[ vol_range ] )
			ax.plot( pred_data[ vol_range ] )

			curr_y_max = np.max( [ np.max( raw_data[ vol_range ] ),
			                       np.max( pred_data[ vol_range ] )
			                     ]
			                   )

			if curr_y_max > y_max:
				y_max = curr_y_max

			curr_y_min = np.min( [ np.min( raw_data[ vol_range ] ),
			                       np.min( pred_data[ vol_range ] )
			                     ]
			                   )

			if curr_y_min < y_min:
				y_min = curr_y_min

			ax_list.append( ax )

		for ( i_run, ax ) in enumerate( ax_list ):

			ax.set_ylim( [ y_min - y_head, y_max + y_head ] )

			_cleanup_fig( ax )

			ax.set_title( "Run %d" % ( i_run + 1 ) )

			ax.set_xlabel( "Time (volumes)" )

			ax.set_ylabel( "Signal (a.u.)" )

		fig.suptitle( roi_name.upper() )

		plt.subplots_adjust( left = 0.09,
		                     bottom = 0.08,
		                     right = 0.98,
		                     top = 0.90,
		                     wspace = 0.39,
		                     hspace = 0.73
		                   )

		if save_path is not None:

			fig_path = "%s-%s.pdf" % ( save_path, roi_name )

			plt.savefig( fig_path )

	if save_path is None:
		plt.show()


def plot_roi_psc( paths, conf ):
	"""Plot the PSC for each ROI"""

	roi_order = [ "V1", "V2", "V3",
	              "V3AB", "LO1", "LO2",
	              "hV4", "VO1", "hMTp"
	            ]

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7.08661, 7.08661, forward = True )

	gs = gridspec.GridSpec( 3, 3 )

	x = np.array( conf[ "stim" ][ "coh_levels" ] ) * 100

	subj_col = [ 0.8 ] * 3

	for ( i_roi, roi_name ) in enumerate( roi_order ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		data = np.loadtxt( "%s-%s.txt" % ( paths[ "roi_mean" ], roi_name.lower() ),
		                   usecols = [ 1, 2, 3, 4 ]
		                 )

		for i_subj in xrange( data.shape[ 0 ] ):

			ax.plot( x,
			         data[ i_subj, : ],
			         color = subj_col
			       )

			ax.scatter( x,
			            data[ i_subj, : ],
			            edgecolor = [ 1 ] * 3,
			            facecolor = subj_col,
			          )

		ax.plot( x,
		         np.mean( data, axis = 0 ),
		         "k",
		         linewidth = 1.5
		       )

		ax.scatter( x,
		            np.mean( data, axis = 0 ),
		            edgecolor = [ 1 ] * 3,
		            facecolor = "k",
		            zorder = 100,
		            marker = "s",
		            s = 35
		       )

		_cleanup_fig( ax )

		ax.set_xlim( [ -10, 110 ] )
		ax.set_ylim( [ -0.4, 0.4 ] )

		ax.set_ylabel( "Response (psc)" )
		ax.set_xlabel( "Stimulus coherence (%)" )

		ax.set_xticks( x )

		ax.text( 0.1,
		         0.9,
		         roi_name,
		         transform = ax.transAxes,
		         fontsize = 10 / 1.25
		       )

	plt.subplots_adjust( left = 0.09,
	                     bottom = 0.07,
	                     right = 0.97,
	                     top = 0.97,
	                     wspace = 0.41,
	                     hspace = 0.34
	                   )

	plt.show()


def plot_lin_trend_hist( paths, conf ):
	"""a"""

	roi_order = [ "V1", "V2", "V3",
	              "V3AB", "LO1", "LO2",
	              "hV4", "VO1", "hMTp"
	            ]

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7.08661, 7.08661, forward = True )

	gs = gridspec.GridSpec( 3, 3 )

	x = np.linspace( -10.0, 10.0, 500 )
	px = 3

	ix = np.logical_and( x > -px, x < px )

	for ( i_roi, roi_name ) in enumerate( roi_order ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		mu = []

		for ( i_subj, subj_id ) in enumerate( conf[ "all_subj" ] ):

			subj_conf = glass_coherence_block.config.get_conf( subj_id )
			subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

			psc = np.vstack( ( np.loadtxt( "%s_%s_%s.txt" % ( subj_paths[ "rois" ][ "psc" ],
			                                                  roi_name.lower(),
			                                                  hemi
			                                                )
			                             )
			                   for hemi in [ "lh", "rh" ]
			                 )
			               )

			lin_coeff = np.sum( psc * np.array( [ -3, -1, +1, +3 ] ), axis = 1 )

			kde = scipy.stats.gaussian_kde( lin_coeff )

			y = kde.evaluate( x )

			y /= np.sum( y )

			ax.plot( x[ ix ], y[ ix ], color = [ 0.25 ] * 3 )

			mu.append( np.mean( lin_coeff ) )

		ylim = ax.get_ylim()

		ax.plot( [ 0, 0 ], [ 0, 0.05 ], "k--", zorder = -100 )

		for subj_mu in mu:

			mu_h = 0.0025

			ax.plot( [ subj_mu ] * 2,
			         [ -mu_h - 0.00125, -mu_h + 0.00125 ],
			         color = [ 0.25 ] * 3
			       )

		_cleanup_fig( ax )

		ax.set_ylim( [ -0.005, 0.05 ] ) # ylim[ 1 ] ] )

		ax.set_xlabel( "Linear trend coefficient" )
		ax.set_ylabel( "Density (norm)" )

		ax.set_yticks( [ 0, 0.01, 0.02, 0.03, 0.04, 0.05 ] )

		ax.text( 0.1,
		         0.9,
		         roi_name,
		         transform = ax.transAxes,
		         fontsize = 10 / 1.25
		       )

	plt.subplots_adjust( left = 0.09,
	                     bottom = 0.07,
	                     right = 0.97,
	                     top = 0.97,
	                     wspace = 0.41,
	                     hspace = 0.34
	                   )

	plt.show()


def plot_trends( paths, conf ):
	"""a"""

	roi_order = [ "V1", "V2", "V3",
	              "V3AB", "LO1", "LO2",
	              "hV4", "VO1", "hMTp"
	            ]

	trends = [ "Linear", "Quadratic", "Cubic" ]

	trend_cf = [ [ -3, -1, +1, +3 ],
	             [ +1, -1, -1, +1 ],
	             [ -1, +3, -3, +1 ]
	           ]

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 3, 7, forward = True )

	gs = gridspec.GridSpec( 3, 1 )

	x_off = 0.2

	for ( i_trend, trend_name ) in enumerate( trends ):

		ax = plt.subplot( gs[ i_trend ] )
		ax.hold( True )

		for ( i_roi, roi_name ) in enumerate( roi_order ):

			stat = np.loadtxt( "%s-%s.txt" % ( paths[ "roi_stat" ], roi_name.lower() ) )
			grp = np.loadtxt( "%s-%s.txt" % ( paths[ "roi_mean" ], roi_name.lower() ),
			                  usecols = np.arange( 1, 5 )
			                )

			for i_subj in xrange( grp.shape[ 0 ] ):
				subj_cf = np.sum( grp[ i_subj, : ] * trend_cf[ i_trend ] )

				ax.scatter( i_roi, subj_cf, color = [ 0.8 ] * 3 )

			ax.plot( [ i_roi - x_off ] * 2,
			         [ stat[ 3, i_trend ], stat[ 4, i_trend ] ],
			         "w"
			       )

			ax.plot( [ i_roi + x_off ] * 2,
			         [ stat[ 3, i_trend ], stat[ 4, i_trend ] ],
			         "w"
			       )

			ax.scatter( i_roi, stat[ 0, i_trend ] )

	plt.show()
