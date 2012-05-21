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

    s1000 = { "subj_id" : "s1000",
              "acq_date" : "20120427",
              "n_runs" : 10,
              "n_fmaps" : 1,
              "run_mot_order" : ( 6, 7, 8, 9, 10, 1, 2, 3, 4, 5 ),
              "x_range" : ( 95, 70 ),
              "y_range" : ( 6, 20 ),
              "z_range" : ( 3, 6 ),
              "comments" : ""
            }

``N_range`` refers to the field-of-view (in voxels) to extract from the images. Use ``fslview`` to check these.


Pre-processing
--------------

Most of the pre-processing is done with the command ``radial_bias_layers_preproc``.
For help on using this script, run::

    radial_bias_layers_preproc --help

Typical usage is::

    radial_bias_layers_preproc sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the preprocessing stage (see below).

The stages are as follows:

Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files and crops the FOV::

    radial_bias_layers_preproc sXXXX convert

After execution, open up each NIFTI file and inspect for image quality.


Correction
~~~~~~~~~~

Applies a motion correction procedure and creates a session mean image::

    radial_bias_layers_preproc sXXXX correct

*N.B. This stage takes quite a while...*

After execution, open up the session summary image that it creates and view in movie mode. This gives a good sense for how well the motion correction worked. You can also inspect the saved motion correction estimates to see how much movement there was.


ROI to images
~~~~~~~~~~~~~

Converts the raw ROI files from mrLoadRet into NIFTI masks::

    radial_bias_layers_preproc sXXXX roi-img

To check this has worked correctly, load the subject's anatomical image and overlay the ROI images - they should lie within expected locations.


Coregistration
~~~~~~~~~~~~~~

The anatomical and ROI images are in a completely different space to the functionals, so they need to be coregistered.

This is somewhat tricky for the high-resolutions we have for the functionals, so it needs to be a multi-step process.

Rough alignment
^^^^^^^^^^^^^^^

The coregistration algorithm is helped enormously if the images are in rough world-space alignment before it begins.

#. In SPM, click ``Display`` and select the **magnitude** fieldmap image.
#. Place the crosshairs over a prominent landmark, such as the furthest posterior region of the occipital lobes. Note down the 3 values in the ``mm`` box.
#. Click ``Display`` again, this time selecting the anatomical image.
#. Place the crosshairs over the same landmark as was used in the magnitude image, and again note the 3 values in the ``mm`` box.
#. Subtract (element-wise) the anatomical ``mm`` values from the magnitude image ``mm`` values, and use the output to populate the ``right``, ``forward``, and ``up`` fields.
#. To check your calculations, change the ``mm`` field to match what it was for the magnitude image and the crosshairs should move to the same landmark.
#. Click ''Reorient images'' and select the anatomical **and the ROI, distance, and ret images**.


Initial coregistration
^^^^^^^^^^^^^^^^^^^^^^

#. In SPM, click ``Coregister (Estimate)``.
#. As the ``Reference image``, select the magnitude image.
#. As the ``Source image``, select the anatomical image.
#. As the ``Other images``, select all the ROI, distance, and ret images.
#. Under ``File``, click ``Save batch`` and call it ``coreg_a.mat`` under the ``anat`` directory.
#. Click on the play icon to set it running.


Final coregistration
^^^^^^^^^^^^^^^^^^^^

#. In SPM, click ``Coregister (Estimate & Reslice)``.
#. As the ``Reference image``, select the mean functional image.
#. As the ``Images to reslice``, select the anatomical image.
#. As the ``Other images``, select the ROI, distance, and ret images.
#. Under ``Estimation options``, change ``Separation`` to ``[2,1]``.
#. Under ``Reslice options``, change ``Interpolation`` to ``Nearest neighbour`` and ``Filename prefix`` to ``rs_``.
#. Under ``File``, click ``Save batch`` and call it ``coreg_b.mat`` under the ``anat`` directory.
#. Click on the play icon to set it running.


Verification
^^^^^^^^^^^^

To check that the coregistration has performed well:

#. In SPM, click ``Check reg``.
#. Select the mean functional image first, and then the (unresliced) anatomical image.
#. Click around some prominent landmarks and check that the two images are in register.


ROI preparation
~~~~~~~~~~~~~~~

Converts the ROI image masks to a set of coordinates, save in numpy format::

    radial_bias_layers_preproc sXXXX roi


Voxel timecourse extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extracts voxel timecourses for each voxel in each ROI::

    radial_bias_layers_preproc sXXXX vtc


Cortical depth extraction
~~~~~~~~~~~~~~~~~~~~~~~~~

Extracts the gray matter depth for each voxel in each ROI::

    radial_bias_layers_preproc sXXXX depth


Retinotopy phase extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extracts the retinotopic phase (corresponding to visual field location) for each voxel in each ROI::

    radial_bias_layers_preproc sXXXX ret


Design
~~~~~~

Computes the experimental design from the logfiles::

    radial_bias_layers_preproc sXXXX design

The extracted design corresponds to the trimmed and HRF corrected voxel timecourses.


Filtering
~~~~~~~~~

Applies a high-pass filter to each voxel in each ROI::

    radial_bias_layers_preproc sXXXX vtc-filt


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
