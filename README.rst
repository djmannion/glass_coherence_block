.. highlight:: bash

=========================================================
Analysis code for Glass coherence block design experiment
=========================================================

Requirements
============

- Python >= 2.7
- numpy
- matplotlib
- scipy
- fmri_tools (`http://www.bitbucket.org/djmannion/fmri_tools <http://www.bitbucket.org/djmannion/fmri_tools/>`_)
- FSL >= 4.1.9
- AFNI/SUMA

Processing stages
=================

Prepare the filesystem
----------------------

1. Make the subject's directory structure::

    mkdir -p sXXXX/{analysis,fmap,func/run{01,02,03,04,05,06,07,08,09,10,11,12},logs,reg,rois}

2. Copy the subject's runtime logfiles to the ``logs`` directory.

3. Make symlinks named ``raw`` in each functional run directory that link to the location of its associated raw DICOM directory::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-such_and_such raw

4. Similarly, make symlinks named ``mag-raw`` and ``ph-raw`` in each fieldmap directory that link to the locations of the fieldmap acquisition::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-SEyada mag-raw
    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/PH-SEyada ph-raw


Update the experiment information file
--------------------------------------

Edit ``get_subj_conf`` within ``glass_coherence_block/config.py`` and add the new subject's information.


Pre-processing
--------------

Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    glass_coherence_block_preproc sXXXX convert

After execution, open up each NIFTI file and inspect for image quality and inspect the summary image to see how much motion there was and as a comparison for the next step.


Correction
~~~~~~~~~~

Applies a motion correction procedure and creates a summary image::

    glass_coherence_block_preproc sXXXX correct

After execution, open up the session summary image that it creates and view in movie mode. This gives a good sense for how well the motion correction worked. You can also inspect the saved motion correction estimates to see how much movement there was.


Fieldmaps
~~~~~~~~~

Prepares the fieldmap::

    glass_coherence_block_preproc sXXXX fieldmap


Unwarping
~~~~~~~~~

Use the fieldmaps to unwarp the functional images to remove the spatial distortion::

    glass_coherence_block_proc sXXXX undistort

To verify that the unwarping has worked correctly:

* Run ``fslview``.
* Load the original or corrected image from a given run.
* Add the magnitude image from the fieldmap as an overlay.
* Notice the geometric distortions in the functional data.
* Add the undistorted image as an overlay, and hide the uncorrected image.
* Toggle the visibility of the undistorted image, and verify that the geometry now aligns well with that of the fieldmap's magnitude image.

Also, look at the session summary image produced and make sure that all looks good across the session.


Coregistration
~~~~~~~~~~~~~~

Make copies of the reference anatomical and mean functional in the ``reg`` directory::

  3dcopy /labs/olmanlab/FsAnatomy/sXXXX/SUMA/sXXXX_SurfVol+orig sXXXX_glass_coherence_block-anat_ref+orig
  3dcopy ../func/sXXXX_glass_coherence_block-mean.nii sXXXX_glass_coherence_block-mean+orig

Run::

    glass_coherence_block_preproc sXXXX sess_reg

Then, with the mean functional as underlay, change the overlay in AFNI to be ``sXXXX_glass_coherence-anat_ref_al+orig`` and confirm that the alignment is perfect.

Finally, set the underlay to be the mean functional and then fire up SUMA::

    suma -spec /labs/olmanlab/FsAnatomy/sXXXX/SUMA/sXXXX_both.spec -sv sXXXX_glass_coherence_block-anat_ref_al+orig

Press ``t`` inside SUMA to send the surfaces to AFNI, and check that they are well registered to the mean functional.

If the alignment is **not** good:

* Fire up ``afni`` from within the ``reg`` directory.
* Set the mean functional as the underlay, and note down the position (in mm) of a landmark in posterior cortex (and the direction indicators).
* Then, set the underlay to the (unaligned) anatomical, and note the position of the same landmark.
* Subtract the two vectors elementwise to give a rough estimate of the translation required, and pass as extra parameters in the subject's configuration.


Surface projection
~~~~~~~~~~~~~~~~~~

The functional images, in their volume space, are now projected onto the cortical surface by averaging between the white matter (smoothed) and pial surfaces::

    glass_coherence_block_preproc sXXXX vol_to_surf


Design preparation
~~~~~~~~~~~~~~~~~~

We need to extract the stimulus and experiment design information from the log files and output it in a format suitable for reading into AFNI's GLM analysis programs::

    glass_coherence_block_preproc sXXXX design_prep


Subject-level analysis
----------------------

GLM
~~~

Runs a GLM analysis::

    glass_coherence_block_subj_analysis sXXXX glm


Localiser mask
~~~~~~~~~~~~~~

Creates a mask of activated nodes from the GLM analysis::

    glass_coherence_block_subj_analysis sXXXX loc_mask


Percent signal change
~~~~~~~~~~~~~~~~~~~~~

Converts the fitted beta values from the GLM to percent signal change::

    glass_coherence_block_subj_analysis sXXXX beta_to_psc


ROI preparation
~~~~~~~~~~~~~~~

First, make copies, within the ``rois`` directory, of the visual localiser ROI datasets::

    3dcopy /labs/olmanlab/FsAnatomy/sXXXX/rois/sXXXX_vis_loc-rois_lh-full.niml.dset sXXXX_glass_coherence_block-vis_loc_rois_lh-full.niml.dset
    3dcopy /labs/olmanlab/FsAnatomy/sXXXX/rois/sXXXX_vis_loc-rois_rh-full.niml.dset sXXXX_glass_coherence_block-vis_loc_rois_rh-full.niml.dset

Guided by the localiser mask and the visual localisers, draw ROIs for the dorsal (value of 100) and ventral (value of 200) 'responsive areas'.
Save these as ``lh_dra``, ``lh_vra``, ``rh_dra``, and ``rh_vra`` in the ``rois`` directory.

Then, run::

    glass_coherence_block_subj_analysis sXXXX roi_prep

Load the resulting ROI dataset (``sXXXX_glass_coherence_block-rois_HEMI-full.niml.dset``) and check that the ROIs are correct (particularly that dorsal and ventral are accurate).


ROI statistics
~~~~~~~~~~~~~~

Extracts the node percent signal change values for each ROI, combines across hemispheres, and computes the contrast coefficents::

    glass_coherence_block_subj_analysis sXXXX roi_xtr


Visual field
~~~~~~~~~~~~

First, draw ROIs for dorsal (101) and ventral (201) cortex, saved as ``lh_dorsal``, ``lh_ventral``, ``rh_dorsal``, ``rh_ventral`` in the ``rois`` directory.
Draw the dorsal ROI first, since it has precedence in overlap.
Then, make sure ``wedge_date`` is set correctly in the subject's configuration.

Then, run::

    glass_coherence_block_subj_analysis sXXXX roi_vf

Task
~~~~

Analyses performance on the behavioural task::

    glass_coherence_block_subj_analysis sXXXX task


Group processing
----------------

ROI preparation
~~~~~~~~~~~~~~~

Aggregates the ROI data for each subject and performs normalisation::

    glass_coherence_block_group_analysis roi_prep


Permutation testing
~~~~~~~~~~~~~~~~~~~

Computes the trend coefficients from the observed data and also generates permutation distributions::

    glass_coherence_block_group_analysis roi_perm


Statistics
~~~~~~~~~~

Generates descriptive and inferential statistics::

    glass_coherence_block_group_analysis roi_stat


Visual field
~~~~~~~~~~~~

As above, but for the effect of visual field position in V3::

    glass_coherence_block_group_analysis vf_v3_prep
    glass_coherence_block_group_analysis vf_v3_perm
    glass_coherence_block_group_analysis vf_v3_stat


Task
~~~~

Runs statistical tests on the task performance across subjects::

    glass_coherence_block_group_analysis task



