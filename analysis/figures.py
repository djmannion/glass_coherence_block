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


def plot_task_perf( conf, paths, show_plot = False ):
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

	if show_plot:
		fig.show()
	else:
		fig_path = paths.fig_task.full( ".svg" )
		plt.savefig( fig_path )



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


def plot_psc( conf, paths, show_plot = False ):
	"""Plot the PSC for each ROI"""

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7.08661, 4.5, forward = True )

	gs = gridspec.GridSpec( 2, 3 )

	x = np.array( conf[ "stim" ][ "coh_levels" ] ) * 100

	subj_col = [ 0.8 ] * 3

	for ( i_roi, ( roi_name, _ ) ) in enumerate( conf[ "ana" ][ "rois" ] ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		subj_data_path = paths.psc.full( "_{roi:s}-norm.txt".format( roi = roi_name ) )

		# this is subjects x coherences
		subj_data = np.loadtxt( subj_data_path )

		for i_subj in xrange( subj_data.shape[ 0 ] ):

			ax.plot( x,
			         subj_data[ i_subj, : ],
			         color = subj_col
			       )

			ax.scatter( x,
			            subj_data[ i_subj, : ],
			            edgecolor = [ 1 ] * 3,
			            facecolor = subj_col,
			          )

		data_path = paths.descrip.full( "_{roi:s}.txt".format( roi = roi_name ) )

		# ( mean, sem )
		( data_mean, data_sem ) = np.loadtxt( data_path )

		ax.plot( x,
		         data_mean,
		         "k",
		         linewidth = 1.5
		       )

		_ = [ ax.plot( [ xx ] * 2,
		               [ xx_data_mean - xx_data_sem, xx_data_mean + xx_data_sem ],
		               "k",
		               linewidth = 1.5
		             )
		      for ( xx, xx_data_mean, xx_data_sem ) in zip( x, data_mean, data_sem )
		    ]

		ax.scatter( x,
		            data_mean,
		            edgecolor = [ 1 ] * 3,
		            facecolor = "k",
		            zorder = 100,
		            marker = "s",
		            s = 35
		       )

		_cleanup_fig( ax )

		ax.set_xlim( [ -10, 110 ] )
		ax.set_ylim( [ -0.325, 0.325 ] )

		if i_roi == 3:
			ax.set_ylabel( "Response (psc)" )
			ax.set_xlabel( "Stimulus coherence (%)" )

		ax.set_xticks( x )
		ax.set_yticks( [ -0.2, 0, 0.2 ] )

		ax.text( 0.1,
		         0.9,
		         conf[ "ana" ][ "roi_labels" ][ i_roi ],
		         transform = ax.transAxes,
		         fontsize = 10 / 1.25
		       )

	plt.subplots_adjust( left = 0.09,
	                     bottom = 0.12,
	                     right = 0.97,
	                     top = 0.97,
	                     wspace = 0.41,
	                     hspace = 0.34
	                   )

	if show_plot:
		plt.show()
	else:
		save_path = paths.fig_psc.full( ".svg" )
		plt.savefig( save_path )


def plot_lin_kde( conf, path, show_plot = False ):
	"""Linear trend histograms"""

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7.08661, 7.08661, forward = True )

	gs = gridspec.GridSpec( 2, 3 )

	x = np.linspace( -10.0, 10.0, 500 )
	px = 3

	ix = np.logical_and( x > -px, x < px )

	for ( i_roi, ( roi_name, roi_id ) ) in enumerate( conf[ "ana" ][ "rois" ] ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		mu = []

		for ( i_subj, subj_id ) in enumerate( conf[ "all_subj" ] ):

			subj_conf = glass_coherence_block.config.get_conf( subj_id )
			subj_paths = glass_coherence_block.analysis.paths.get_subj_paths( subj_conf )

			coef = np.loadtxt( subj_paths.roi.con_coef.full( ".txt" ) )

			i_coef_roi = ( coef[ :, 0 ] == int( roi_id ) )

			lin_coeff = coef[ i_coef_roi, 1 ]

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
