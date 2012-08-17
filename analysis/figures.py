"""
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats


def _set_defaults():
	"""Set some sane defaults for figures.
	"""

	params = { 'axes.labelsize': 8,
	           'axes.titlesize' : 10,
	           'font.family' : 'Arial',
	           'font.sans-serif' : 'Helvetica',
	           'text.fontsize': 12,
	           'legend.fontsize': 7,
	           'xtick.labelsize': 5,
	           'xtick.direction' : 'out',
	           'xtick.major.size' : 2,
	           'ytick.labelsize': 5,
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
	"""Plot the PSC for each ROI, plus CIs and subjects"""

	roi_order = [ "V1", "V2", "V3",
	              "V3AB", "LO1", "LO2",
	              "hV4", "VO1", "hMTp"
	            ]

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7, 7, forward = True )

	gs = gridspec.GridSpec( 3, 3 )

	x = np.array( conf[ "stim" ][ "coh_levels" ] ) * 100

	for ( i_roi, roi_name ) in enumerate( roi_order ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		data = np.loadtxt( "%s-%s.txt" % ( paths[ "roi_mean" ], roi_name.lower() ),
		                   usecols = [ 1, 2, 3, 4 ]
		                 )

		data = ( data.T - np.mean( data, axis = 1 ) ).T

		ci = np.loadtxt( "%s_%s.txt" % ( paths[ "roi_ci" ], roi_name.lower() ) )

		for i_subj in xrange( data.shape[ 0 ] ):

			ax.plot( x,
			         data[ i_subj, : ],
			         color = [ 0.7, 0.7, 0.7 ]
			       )

		assert( np.allclose( np.mean( data, axis = 0 ), ci[ 0, :4 ] ) )

		ax.plot( x,
		         np.mean( data, axis = 0 ),
		         color = "k",
		         linewidth = 3
		       )

		for ci_y in ci[ 1:, :4 ]:

			ax.plot( x,
			         ci_y,
			         "k--"
			       )

		_cleanup_fig( ax )

		ax.set_xlim( [ -10, 110 ] )
		ax.set_ylim( [ -0.4, 0.4 ] )

		ax.set_ylabel( "Response (psc)" )
		ax.set_xlabel( "Coherence" )

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


def plot_trends_old( paths, conf ):
	"""Plot the trends for each ROI"""

	roi_order = [ "V1", "V2", "V3",
	              "V3AB", "LO1", "LO2",
	              "hV4", "VO1", "hMTp"
	            ]

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7, 7, forward = True )

	gs = gridspec.GridSpec( 3, 3 )

	x = np.arange( 3 )

	for ( i_roi, roi_name ) in enumerate( roi_order ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		ax.plot( [ -0.5, 2.5 ], [ 0, 0 ], "k--" )

		ci = np.loadtxt( "%s_%s.txt" % ( paths[ "roi_ci" ], roi_name.lower() ) )

		for ( i_plt, i_ci ) in enumerate( np.arange( 4, 7 ) ):

			ax.plot( [ i_plt, i_plt ],
			         ci[ 1:, i_ci ],
			         'k'
			       )

			ax.scatter( i_plt,
			            ci[ 0, i_ci ]
			          )


		_cleanup_fig( ax )

		ax.set_xlim( [ -0.5, 2.5 ] )
		ax.set_ylim( [ -2, 2 ] )

		ax.set_ylabel( "Trend coefficient" )
		ax.set_xlabel( "Trend" )

		ax.set_xticks( x )
		ax.set_xticklabels( [ "Linear", "Quad.", "Cubic" ] )

		ax.text( 0.8,
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

def plot_vox_psc( paths, conf ):
	"""a"""

	_set_defaults()

	left, width = 0.1, 0.8
	bottom, height = 0.1, 0.6
	bottom_h = left_h = left+height+0.02

	rect_scatter = [ left, bottom, width, height ]
	rect_hist = [ left, bottom_h, width, 0.23 ]

	i_loc_t = np.array( [ 2, 5 ] )

	i_exp_var = 0

	crit_t = scipy.stats.t.ppf( 0.999, 292 )

	dot_cols = [ [ 0, 0, 1 ], [ 1, 0, 0 ] ]

	for ( i_roi, roi ) in enumerate( conf[ "ana" ][ "rois" ] ):

		( roi_name, _ ) = roi

		fig = plt.figure()
		fig.set_size_inches( 8, 5, forward = False )

		ax_scatter = plt.axes( rect_scatter)
		ax_hist = plt.axes( rect_hist )

		xlims = []
		ylims = []

		axes = []

		loc_psc = []
		exp_psc = []
		exp_var = []
		cols = []

		sig = []

		for hemi in [ "lh", "rh" ]:

			exp_stat_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "exp_roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			exp_roi_stat = np.loadtxt( exp_stat_file )

			loc_stat_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "loc_roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			loc_roi_stat = np.loadtxt( loc_stat_file )

			loc_psc_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "loc_roi_psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			loc_roi_psc_all = np.loadtxt( loc_psc_file )

			i_psc = np.argmax( loc_roi_stat[ :, i_loc_t ], axis = 1 )
			t_max = np.max( loc_roi_stat[ :, i_loc_t ], axis = 1 )

			loc_roi_psc = np.empty( loc_roi_psc_all.shape[ 0 ] )

			for i_vox in xrange( loc_roi_psc_all.shape[ 0 ] ):

				loc_roi_psc[ i_vox ] = loc_roi_psc_all[ i_vox, i_psc[ i_vox ] ]

				if t_max[ i_vox ] >= crit_t:
					cols.append( dot_cols[ 0 ] )
					sig.append( 1 )
				else:
					cols.append( dot_cols[ 1 ] )
					sig.append( 0 )

			loc_psc.append( loc_roi_psc.copy() )

			exp_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "exp_roi_psc" ],
			                              roi_name,
			                              hemi
			                            )

			exp_roi = np.loadtxt( exp_file )

			exp_psc.append( exp_roi.copy() )

			exp_var.append( exp_roi_stat[ :, i_exp_var ].copy() )

		loc_psc = np.concatenate( loc_psc )
		exp_psc = np.concatenate( exp_psc )
		exp_var = np.concatenate( exp_var )

		sig = np.array( sig )

		ns_kde = scipy.stats.gaussian_kde( exp_psc[ sig == 0 ] )
		s_kde = scipy.stats.gaussian_kde( exp_psc[ sig == 1 ] )

		ax_scatter.scatter( exp_psc,
		                    loc_psc,
		                    facecolor = cols,
		                    edgecolor = cols,
#		                    s = 1. / exp_var * 10000,
		                    alpha = 0.1
		                  )

		xl = ax_scatter.get_xlim()

		xl_max = np.max( np.abs( xl ) )

		ax_scatter.set_xlim( ( -xl_max, xl_max ) )

		ax_scatter.set_xlabel( "Coherent scene beta (psc)" )
		ax_scatter.set_ylabel( "Localiser beta (psc)" )

#		ax_hist.hist( exp_psc[ sig == 0 ], orientation = "horizontal" )

		x = np.linspace( -xl_max, xl_max, 500 )

		ns_eval = ns_kde( x )
		s_eval = s_kde( x )

		ax_hist.plot( x, ns_eval, color=dot_cols[ 1 ] )
		ax_hist.hold( True )
		ax_hist.plot( x, s_eval, color=dot_cols[ 0 ] )

		ax_hist_ylim = ax_hist.get_ylim()

		ax_hist.plot( [ 0, 0 ], ax_hist_ylim, "k--" )

		ns_mean_psc = np.mean( exp_psc[ sig == 0 ] )
		s_mean_psc = np.mean( exp_psc[ sig == 1 ] )

		ax_hist.plot( [ ns_mean_psc, ns_mean_psc ],
		              ax_hist_ylim,
		              color = dot_cols[ 1 ],
		              linestyle = "--"
		            )

		ax_hist.plot( [ s_mean_psc, s_mean_psc ],
		              ax_hist_ylim,
		              color = dot_cols[ 0 ],
		              linestyle = "--"
		            )


#			ax.hold( True )

		plt.title( roi_name.upper() )

#		ylims.append( np.max( np.abs( ax.get_ylim() ) ) )

#		axes.append( ax )

#		ylim = np.max( ylims )

#		for ax in axes:

#			xlim = ax.get_xlim()

#			ax.plot( xlim, [ 0, 0 ], "k--" )

#			ax.set_ylim( ( -ylim, ylim  ) )
#			ax.set_xlim( xlim )


#		plt.savefig( "/home/dmannion/im_temp/ns_ap_%s_%s.png" % ( conf[ "subj" ][ "subj_id" ], roi_name ) )
	plt.show()


def subj_cond_diff( paths, conf ):
	"""
	"""

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 4.8, 3, forward = True )

	ax = fig.gca()

	ax.hold( True )

	for ( i_roi, roi_name ) in enumerate( conf[ "ana" ][ "rois" ] ):

		diff = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "cond_diff" ],
		                                roi_name
		                              )
		              )

		ax.plot( ( i_roi, i_roi ), diff[ 1: ], "k" )

		ax.scatter( i_roi, diff[ 0 ],
		            facecolor = "k",
		            edgecolor = "w",
		            s = 35,
		            zorder = 100
		           )

	ax.plot( ax.get_xlim(), ( 0, 0 ), "k--" )

	_cleanup_fig( ax )

	ax.set_xlim( ( -0.5, len( conf[ "ana" ][ "rois" ] ) - 0.5 ) )

	ax.set_ylim( ( -0.5, 0.5 ) )
	ax.set_xticks( np.arange( len( conf[ "ana" ][ "rois" ] ) ) )
	ax.set_xticklabels( conf[ "ana" ][ "rois" ] )

	ax.set_ylabel( "Same - different scene (psc)" )
	ax.set_xlabel( "Visual area" )

	fig = plt.gcf()
	fig.set_size_inches( 4.8, 3 )

	fig.show()

	fig = plt.gcf()
	fig.set_size_inches( 4.8, 3 )

	plt.show()


def subj_resp_over_time( paths, conf ):
	"""
	"""

	_set_defaults()



	fig = plt.figure()

	fig.set_size_inches( 7, 5, forward = True )

	gs = gridspec.GridSpec( 2, 3 )

	for ( i_roi, roi_name ) in enumerate( conf[ "ana" ][ "rois" ] ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		blks = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "block" ],
		                                roi_name
		                              )
		              )

		cond_a = blks[ blks[ :, 1 ] == 0, 0 ]
		cond_b = blks[ blks[ :, 1 ] == 1, 0 ]


		ax.plot( cond_a, "b" )
		ax.plot( cond_b, "r" )

		ax.plot( ax.get_xlim(), ( 0, 0 ), "k--", zorder = -1 )

		ax.text( 0.82,
		         0.9,
		         roi_name,
		         transform = ax.transAxes,
		         fontsize = 10 / 1.25
		       )

		_cleanup_fig( ax )

		ax.set_ylabel( "Response (psc)" )
		ax.set_xlabel( "Time (blocks)" )

	plt.subplots_adjust( left = 0.10,
	                     bottom = 0.10,
	                     right = 0.97,
	                     top = 0.90,
	                     wspace = 0.40,
	                     hspace = 0.34
	                   )

	plt.show()
