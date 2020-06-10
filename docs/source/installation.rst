#####################
Installation
#####################

PPMAP is currently configured to run on the Super Computing Wales Hawk architecture which runs the Slurm scheduler, and the intel ifort compiler. The algorithm is hosted on GitHub:

 `https://github.com/kennethmarsh/ppmap <https://github.com/kennethmarsh/ppmap>`_

=========================
Installing the Algorithm
=========================

To install:

* Clone the repository to the desired directory in your cluster ``/home/[USERNAME]/`` space.
.. note:: After cloning, you should have a ``ppmap/`` directory with the path ``/home/[USERNAME]/ppmap/``, which contains a futher ``source/`` directory and ``templates/`` directory. If you cloned the repository in a location other than  ``/home/[USERNAME]/``, make a note of this address **without** the ``ppmap/`` direcotry (e.g. ``/home/[USERNAME]/some/other/path/``) and see the **Warning** below.
* Open the makefile in a text editor and edit the ``USRNAME`` variable with your Hawk username (e.g. c.c1234567).
.. warning:: If the ``ppmap`` directory is located somewhere other than ``/home/[USERNAME]/``, you will need to go into the ``/ppmap/templates/`` directory and open ``template_run_ppmap.sh`` in a text editor. You will then need to edit the line ``code=`` variable from ``code=${HOME}/ppmap/bin/ppmap`` to include the directory you copied above (e.g. ``code=/home/[USERNAME]/some/other/path/ppmap/bin/ppmap``). Do not be concerned that the ``ppmap/`` directory does not contain a ``bin/`` directory yet. One will be created during the installation process.
* Move into the ``ppmap/`` directory and run the makefile:

.. code-block:: shell

	make
    
This should compile the PreMAP and PPMAP executables and place them in the ``ppmap/bin/`` directory. It should also generate the ``run_[FIELD]_all`` and ``run_ppmap.sh`` script files from the templates, configure them for your username, and place them into the ``ppmap/`` directory. 

PPMAP should now be compiled and ready to operate.