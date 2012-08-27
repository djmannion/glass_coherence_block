=========================================================
Analysis code for Glass coherence block design experiment
=========================================================

Requirements
============

- Python >= 2.7
- numpy
- matplotlib
- scipy
- scikits.bootstrap
- fmri_tools (`http://www.bitbucket.org/djmannion/fmri_tools <http://www.bitbucket.org/djmannion/fmri_tools/>`_)
- FSL >= 4.1.9
- AFNI/SUMA
- SPM

Processing stages
=================

Prepare the filesystem
----------------------

1. Make the subject's directory structure::

    mkdir -p sXXXX/{analysis,fmap/f01,func/run{01,02,03,04,05,06,07,08,09,10,11,12},logs,reg,rois}

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

Most of the pre-processing is done with the command ``glass_coherence_block_proc``.
For help on using this script, run::

    glass_coherence_block_proc --help

Typical usage is::

    glass_coherence_block_proc sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the processing stage (see below).

The stages are as follows:

Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    glass_coherence_block_proc sXXXX convert

After execution, open up each NIFTI file and inspect for image quality and inspect the summary image to see how much motion there was and as a comparison for the next step.


Correction
~~~~~~~~~~

Applies a motion correction procedure and creates session mean and summary images::

    glass_coherene_block_proc sXXXX correct

*N.B. This stage takes quite a while...*

After execution, open up the session summary image that it creates and view in movie mode. This gives a good sense for how well the motion correction worked. You can also inspect the saved motion correction estimates to see how much movement there was.


Fieldmaps
~~~~~~~~~

Prepares the fieldmaps::

    glass_coherence_block_proc sXXXX fieldmap


Unwarping
~~~~~~~~~

Before running, need to have made a symbolic link in each functional run directory to that run's fieldmap. For example::

    ln -s ../../fmap/f1/sXXXX_glass_coherence_block_fmap_1-fmap.nii sXXXX_glass_coherence_block_run_1-fmap.nii

Then, to use the fieldmaps to unwarp the functional images to remove the spatial distortion::

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

The anatomical image is in a completely different space to the functionals, so they need to be coregistered.

The automatic FSL tools are *horrible* at doing this coregistration (in my experience), so we need to do it more manually using SPM.

Follow the instructions `here <http://visual-localiser-analysis-notes.readthedocs.org/en/latest/func.html#coregister-base-anatomy-to-functional-session>`__ to get the registration, then run::

    glass_coherence_block_proc sXXXX surf_reg


Surface projection
~~~~~~~~~~~~~~~~~~

The functional images, in their volume space, are now projected onto the cortical surface by averaging between the white matter and pial surfaces::

    glass_coherence_block_proc sXXXX vol_to_surf


Design preparation
~~~~~~~~~~~~~~~~~~

We need to extract the stimulus and experiment design information from the log files and output it in a format suitable for reading into AFNI's GLM analysis programs::

    glass_coherence_block_proc sXXXX design_prep


Subject-level analysis
----------------------

GLM
~~~

Runs a GLM analysis::

    glass_coherence_block_proc sXXXX glm


Localiser mask
~~~~~~~~~~~~~~

Creates a mask of activated nodes from the GLM analysis::

    glass_coherence_block_proc sXXXX loc_mask


Percent signal change
~~~~~~~~~~~~~~~~~~~~~

Converts the fitted beta values from the GLM to percent signal change::

    glass_coherence_block_proc sXXXX beta_to_psc


ROI statistics
~~~~~~~~~~~~~~

Extracts the node values for each ROI::

    glass_cohernce_block_proc sXXXX roi_xtr


Adjusted timecourses
~~~~~~~~~~~~~~~~~~~~

Synthesizes raw and predicted timecourses that are adjusted to remove baseline trends::

    glass_coherence_block_proc sXXXX raw_adj


ROI timecourses
~~~~~~~~~~~~~~~

Compiles raw and predicted timecourses (adjusted) for each ROI::

    glass_coherence_block_proc sXXXX roi_tc


Timecourse visualisation
~~~~~~~~~~~~~~~~~~~~~~~~

Plots the average raw and predicted timecourses (adjusted) for each run (panel) and ROI (figure)::

    glass_coherence_block_proc sXXXX plot_tc


Group processing
----------------

ROI aggregation
^^^^^^^^^^^^^^^

Grabs the ROI data for each subject/hemisphere and collates it together by averaging over nodes::

  glass_coherence_block_group rois


Statistics
^^^^^^^^^^

Compute the trend coefficients and run permutation tests to get significance values (p and q)::

  glass_coherence_block_group stat











