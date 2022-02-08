=============
jjsoft plugin
=============

This plugin allows to use jjsoft within the Scipion framework. The binaries for the plugin must be obtained from 
3demimageprocessing_.

=====
Setup
=====

- **Install this plugin:**

.. code-block::

    scipion installp -p scipion-em-jjsoft

OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

Alternatively, in devel mode:

.. code-block::

    scipion installp -p local/path/to/scipion-em-jjsoft --devel

===================================
Manual installation of the binaries
===================================

Four binaries are used by jjsoft protocols. All of them can be found in 3demimageprocessing_. These are:

- tomo3d_
- tomowarpalign_
- tomobflow_
- tomoeed_

To get each of them, click on the corresponding link and fill the form located at the bottom of the displayed webpage.
Once all of them have been downloaded, place them in scipionDirectory/software/em/jjsoft-VersionNumber. To avoid
version binary names mismatching, they should be renamed or linked, so the names that appear inside the mentioned
directory are tomo3d, tomowarpalign, tomobflow and tomoeed, respectively. Thus, the final file structure expected by
Scipion should look like this:

scipionDirectory/software/em/jjsoft-VersionNumber
-tomo3d
-tomowarpalign
-tomobflow
-tomoeed
-jjsoft_installed (plugin installation file, generated automatically during the plugin installation)



.. _3demimageprocessing: https://sites.google.com/site/3demimageprocessing/
.. _tomo3d: https://sites.google.com/site/3demimageprocessing/tomo3d
.. _tomowarpalign: https://sites.google.com/site/3demimageprocessing/tomoalign
.. _tomobflow: https://sites.google.com/site/3demimageprocessing/tomobflow
.. _tomoeed: https://sites.google.com/site/3demimageprocessing/tomoeed