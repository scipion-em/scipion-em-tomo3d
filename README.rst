=======================
Scipion - PySeg plugin
=======================

This plugin allows to use some of the software packages from 3demimageprocessing_ within the Scipion framework.
**The binaries for the plugin MUST be downloaded manually**:

====================================
Setup - Binaries manual installation
====================================

The following binaries are used by tomo3d protocols. All of them can be found in 3demimageprocessing_. These are:

- tomo3d_
- tomobflow_
- tomoeed_

To get each of them, click on the corresponding link and fill the form located at the bottom of the displayed webpage.
Once all of them have been downloaded, place them in scipionDirectory/software/em/tomo3d-VersionNumber.Thus, the final 
file structure expected by Scipion should look like this, being the listed files below the binaries or a link pointing to each:

scipionDirectory/software/em/tomo3d-VersionNumber

- tomo3d
- tomobflow
- tomoeed

============
Installation
============

The plugin can be installed in user (stable) or developer (latest, may be unstable) mode:

**1. User (stable) version:**:

.. code-block::

    scipion installp -p scipion-em-tomo3d

**2. Developer (latest, may be unstable) version:**:

* Clone the source code repository:

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-tomo3d.git
    
* Install:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-tomo3d --devel

=========
Protocols
=========
The integrated protocols are:

1. tomo3d - denoise tomogram:

2. tomo3d - reconstruct tomogram:

=====
Tests
=====

The installation can be checked out running some tests:

.. code-block::

     scipion3 tests tomo3d.tests.test_protocol_reconstruct_tomogram

.. code-block::

    scipion3 tests tomo3d.tests.test_protocol_denoise_tomogram
    

==========
References
==========

* `Fast tomographic reconstruction on multicore computers. <https://doi.org/10.1093/bioinformatics/btq692>`_
  J.I. Agulleiro, J.J. Fernandez. Bioinformatics 27:582-583, 2011.

* `Tomo3D 2.0 – Exploitation of Advanced Vector eXtensions (AVX) for 3D reconstruction. <https://doi.org/10.1016/j.jsb.2014.11.009>`_
  J.I. Agulleiro, J.J. Fernandez. Journal of Structural Biology 189:147-152, 2015.
  
* `TOMOBFLOW: feature-preserving noise filtering for electron tomography. <https://doi.org/10.1186/1471-2105-10-178>`_
  J.J. Fernandez. BMC Bioinformatics 2009, 10:178.
  
* `TomoEED: Fast Edge-Enhancing Denoising of Tomographic Volumes. <https://doi.org/10.1093/bioinformatics/bty435>`_
  J.J. Moreno, A. Martinez-Sanchez, J.A. Martinez, E.M. Garzon, J.J. Fernandez. Bioinformatics 34:3776-3778, 2018. 
  
===================
Contact information
===================

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net, scipion@cnb.csic.es or open an issue_.

We'll be pleased to help.

*Scipion Team*
  

.. _3demimageprocessing: https://sites.google.com/site/3demimageprocessing/
.. _tomo3d: https://sites.google.com/site/3demimageprocessing/tomo3d
.. _tomowarpalign: https://sites.google.com/site/3demimageprocessing/tomoalign
.. _tomobflow: https://sites.google.com/site/3demimageprocessing/tomobflow
.. _tomoeed: https://sites.google.com/site/3demimageprocessing/tomoeed
.. _issue: https://github.com/scipion-em/scipion-em-tomo3d/issues
