#####################
Installation
#####################

PPMAP is currently configured to run on the `Super Computing Wales Hawk Supercomputer <https://portal.supercomputing.wales/index.php/about-hawk/>`_ cluster. For details on the scheduler and compiler, see :ref:`Requirements`. The algorithm is hosted on GitHub:

 `https://github.com/kennethmarsh/ppmap <https://github.com/kennethmarsh/ppmap>`_
 
============
Requirements
============

The current implementation of PPMAP is designed to work with the SLURM scheduler. PPMAP takes advantage of the current Hawk ``HTC`` nodes, which each contain 40 Intel Xeon Skylake Gold 6148 2.4GHz cores and 4.8GB memory per core (for a total of 192GB per node). When PPMAP is run, it breaks down the required task into up to 10 separate jobs, scheduling each on an independent node to run simultaniously. Individual jobs are independent of one another, and thus do not require Open MPI. The algorithm uses all 40 cores on a given node, utilising OpenMP to decrease the operation time of individual submitted jobs. Therefore PPMAP requires OpenMP to function.

PPMAP is currently compiled with the Intel ``ifort`` compiler, and utilises ``-O3`` optimisation. PPMAP is known to perform well with ``ifort Version 10.0.2.199 Build 20180210``.

=========================
Installing the Algorithm
=========================

To install:

* Clone the repository to the desired directory in your cluster ``/home/[USERNAME]/`` space.
.. note:: After cloning, you should have a ``ppmap/`` directory with the path ``/home/[USERNAME]/ppmap/``, which contains a further ``source/`` directory and ``templates/`` directory. If you cloned the repository in a location other than  ``/home/[USERNAME]/``, make a note of this address **without** the ``ppmap/`` direcotry (e.g. ``/home/[USERNAME]/some/other/path/``) and see the **Warning** below.
* Open the makefile in a text editor and edit the ``USRNAME`` variable with your Hawk username (e.g. c.c1234567).
.. note:: If necessary, also edit the ``LOAD`` command to reflect loading the intel compiler on your cluster. If you edit this line, also edit the corrisponding lines in the ``template_run_[FIELD]_all`` file in the ``/ppmap/templates/`` directory.
.. warning:: If the ``ppmap`` directory is located somewhere other than ``/home/[USERNAME]/``, you will need to go into the ``/ppmap/templates/`` directory and open ``template_run_ppmap.sh`` in a text editor. You will then need to edit the line ``code=`` variable from ``code=${HOME}/ppmap/bin/ppmap`` to include the directory you copied above (e.g. ``code=/home/[USERNAME]/some/other/path/ppmap/bin/ppmap``). Do not be concerned that the ``ppmap/`` directory does not contain a ``bin/`` directory yet. One will be created during the installation process.
* Move into the ``ppmap/`` directory and run the makefile:

.. code-block:: console

	$ make
    
This should compile the PreMAP and PPMAP executables and place them in the ``ppmap/bin/`` directory. It should also generate the ``run_[FIELD]_all`` and ``run_ppmap.sh`` script files from the templates, configure them for your username, and place them into the ``ppmap/`` directory. 

PPMAP should now be compiled and ready to operate.