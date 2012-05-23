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

Processing stages
=================

Prepare the filesystem
----------------------

1. Make the subject's directory structure::

    mkdir -p sXXXX/{anat,analysis,fmap/f1,func/run{01,02,03,04,05,06,07,08,09,10,11,12},log,roi}

2. Copy the subject's runtime logfiles to the ``log`` directory.

3. Make symlinks named ``raw`` in each functional run directory that link to the location of its associated raw DICOM directory::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-such_and_such raw

4. Similarly, make symlinks named ``mag-raw`` and ``ph-raw`` in each fieldmap directory that link to the locations of the fieldmap acquisition::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-SEyada mag-raw
    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/PH-SEyada ph-raw

5. Copy (not symlink) the subject's main high-res anatomical (skull stripped) from the main repository to the ``anat`` directorym using FSL::

    fslmaths /labs/olmanlab/Anatomy/sXXXX/sXXXX_stripped sXXXX_anat

  N.B. We want to use single NIFTI files, so before running the above you may need to run::

      setenv FSLOUTPUTTYPE NIFTI

6. Copy the relevant ROI MAT files from the visual localisers repository (the Gray view) to the ``roi`` directory.


Update the experiment information file
--------------------------------------

Edit ``get_subj_conf`` within ``glass_coherence_block/config.py`` and add the new subject's information.

For example::

    s1000 = { "subj_id" : "s1021",
              "acq_date" : "20120522",
              "n_runs" : 12,
              "n_fmaps" : 1,
              "run_st_mot_order" : ( 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6 ),
              "comments" : ""
            }


Pre-processing
--------------

Most of the pre-processing is done with the command ``glass_coherence_block_preproc``.
For help on using this script, run::

    glass_coherence_block_preproc --help

Typical usage is::

    glass_coherence_block_preproc sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the preprocessing stage (see below).

The stages are as follows:

Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    glass_coherence_block_preproc sXXXX convert

After execution, open up each NIFTI file and inspect for image quality and inspect the summary image to see how much motion there was and as a comparison for the next step.


Correction
~~~~~~~~~~

Applies a motion correction procedure and creates session mean and summary images::

    glass_coherene_block_preproc sXXXX correct

*N.B. This stage takes quite a while...*

After execution, open up the session summary image that it creates and view in movie mode. This gives a good sense for how well the motion correction worked. You can also inspect the saved motion correction estimates to see how much movement there was.


Fieldmaps
~~~~~~~~

Prepares the fieldmaps::

    glass_coherence_block_preproc SXXXX fieldmaps


Unwarping
~~~~~~~~

Before running, need to have made a symbolic link in each functional run directory to that run's fieldmap. For example::

    ln -s ../../fmap/f1/sXXXX_glass_coherence_block_fmap_1-fmap.nii sXXXX_glass_coherence_block_run_1-fmap.nii

Then, to use the fieldmaps to unwarp the functional images to remove the spatial distortion::

    glass_coherence_block_preproc sXXXX undistort

To verify that the unwarping has worked correctly:

* Run ``fslview``.
* Load the original or corrected image from a given run.
* Add the magnitude image from the fieldmap as an overlay.
* Notice the geometric distortions in the functional data.
* Add the undistorted image as an overlay, and hide the uncorrected image.
* Toggle the visibility of the undistorted image, and verify that the geometry now aligns well with that of the fieldmap's magnitude image.

Also, look at the session summary image produced and make sure that all looks good across the session.


ROI to images
~~~~~~~~~~~~~

Converts the raw ROI files from mrLoadRet into NIFTI masks::

    glass_coherence_block_preproc sXXXX roi-img

To check this has worked correctly, load the subject's anatomical image and overlay the ROI images - they should lie within expected locations.


Coregistration
~~~~~~~~~~~~~~

The anatomical and ROI images are in a completely different space to the functionals, so they need to be coregistered.

The automatic FSL tools are *horrible* at doing this coregistration (in my experience), so we need to do it more manually using SPM.

Rough alignment
^^^^^^^^^^^^^^^

The coregistration algorithm is helped enormously if the images are in rough world-space alignment before it begins.

#. In SPM, click ``Display`` and select the mean functional image.
#. Place the crosshairs over a prominent landmark, such as the furthest posterior region of the occipital lobes. Note down the 3 values in the ``mm`` box.
#. Click ``Display`` again, this time selecting the anatomical image.
#. Place the crosshairs over the same landmark as was used in the functionals, and again note the 3 values in the ``mm`` box.
#. Subtract (element-wise) the anatomical ``mm`` values from the functional ``mm`` values, and use the output to populate the ``right``, ``forward``, and ``up`` fields.
#. To check your calculations, change the ``mm`` field to match what it was for the functional and the crosshairs should move to the same landmark.
#. Click ''Reorient images'' and select the anatomical **and the ROI images**.

Coregistration
^^^^^^^^^^^^^^

#. In SPM, click ``Coregister (Estimate & Reslice)``.
#. As the ``Reference image``, select the mean functional image.
#. As the ``Images to reslice``, select the anatomical image.
#. As the ``Other images``, select all the ROI images.
#. Under ``Reslice options``, change ``Interpolation`` to ``Nearest neighbour`` and ``Filename prefix`` to ``rs``.
#. Under ``File``, click ``Save batch`` and call it ``coreg.mat`` under the ``anat`` directory.
#. Click on the play icon to set it running.

Verification
^^^^^^^^^^^^

To check that the coregistration has performed well:

#. In SPM, click ``Check reg``.
#. Select the mean functional image first, and then the (unresliced) anatomical image.
#. Click around some prominent landmarks and check that the two images are in register.


ROI preparation
~~~~~~~~~~~~~~

Converts the ROI image masks to a set of coordinates, save in numpy format::

    glass_coherence_block_preproc sXXXX roi


Voxel timecourse extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~

Extracts voxel timecourses for each voxel in each ROI::

    glass_coherence_block_preproc sXXXX vtc

The resulting timecourses have been trimmed.


Voxel culling
~~~~~~~~~~~~~

Removes voxels that have high mean-normalised variance::

    glass_coherence_block_preproc sXXXX vtc-cull


Timecourse averaging
~~~~~~~~~~~~~~~~~~~~

Averages over the voxels in each ROI::

    glass_coherence_block_preproc sXXXX vtc-avg


Design
~~~~~

Computes the experimental design from the logfiles::

    glass_coherence_block_preproc sXXXX design

The extracted design corresponds to the trimmed voxel timecourses.




Subject-level analysis
----------------------

The subject-level analysis is done with the command ``radial_bias_layers_subj_analysis``
For help on using this script, run::

    radial_bias_layers_subj_analysis --help

Typical usage is::

    radial_bias_layers_subj_analysis sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the preprocessing stage (see below).

The stages are as follows:

Blocks
~~~~~~

Extracts the block responses for each condition and ROI::

    radial_bias_layers_subj_analysis sXXXX blocks



Analysis datafiles
==================

The pre-processing / analysis pipeline produces the following files:

coords-ROI
  ( 3 axes, n voxels ) array of coordinate locations.

coords_sel-ROI
  ( 3 axes, n(s) voxels ) array of coordinate locations, *after* voxel selection based on the localiser analysis.

vtc-ROI
  ( 128 volumes, 10 runs, n voxels ) array of BOLD signals. These are in scanner units, in a timeseries that has been trimmed and HRF corrected.

vtc_sel-ROI
  ( 128 volumes, 10 runs, n(s) voxels ) array of BOLD signals. As above, but only including selected voxels.

loc_vtc_sel-ROI
  ( 128 volumes, 2 runs, n(s) voxels ) array of BOLD signals. As above, but for the localiser data.

vtc_avg-ROI
  ( 128 volumes, 10 runs ) array of BOLD signals. ROI timecourses averaged across all *selected* voxels, high-pass filtered, and covert to percent signal change.

loc_vtc-ROI
  ( 128 volumes, 2 runs, n voxels ) array of BOLD signals. As above, but for the localiser data.

loc_stat-ROI
  ( n voxels, [ t statistic, p value ] ) array of statistics data. These report the results of a left side stimulation > right side stimulation localiser analysis.

design
  ( 16 blocks, 10 runs, [ i_vol, i_cond ) integer array.
  ``i_vol`` is the volume index for the start of the block in a timecourse that has been trimmed and HRF corrected, and ``i_cond`` is the condition.

loc_design
  ( 16 blocks, 2 runs, [ i_vol, i_cond ] ) integer array.
  As above, but for the localiser data.

block
  ( 160 blocks, [ psc, cond, block in run, run ] ) array. Shows the percent signal change of each block, obtained by averaging all the timepoints corresponding to the block.
