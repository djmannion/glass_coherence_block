"""
Set of routines to analyse single-subject fMRI data for the Glass coherence
block design fMRI experiment.
"""

from __future__ import division

import numpy as np

import fmri_tools.glm

def get_glm_beta( paths, conf ):
	"""Fits a GLM and extracts the beta weights for each condition"""

	evt_info = np.load( paths[ "design" ][ "evt_info" ] )

	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the average vtc for this roi
		vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_avg" ],
		                               roi_name
		                             )
		             )

		glm = fmri_tools.glm.RegressGLM( evt_info,
		                                 vtc,
		                                 conf[ "acq" ][ "tr_s" ],
		                                 conf[ "ana" ][ "poly_ord" ]
		                               )

		glm.fit()

		glm_beta = glm.get_beta()

		glm_beta = np.mean( glm_beta, axis = 1 )

		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "beta" ], roi_name ),
		         glm_beta
		       )
