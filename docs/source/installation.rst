#####################
Installation
#####################

PPMAP is currently configured to run on the Super Computing Wales Hawk architecture which runs the Slurm scheduler, and the intel ifort compiler. The algorithm is hosted on GitHub:

 `https://github.com/kennethmarsh/ppmap <https://github.com/kennethmarsh/ppmap>`_

=========================
Compiling the algorithm
=========================

To install, clone the repository to the desired directory and run the makefile

.. code-block:: shell

    make
    
This will compile the PreMAP and PPMAP executables. 

===============================
Configuring the shell scripts
===============================

* Open the ``run_ppmap.sh`` file in a text editor 
* Replace all instances of ``[USERNAME]`` with your Hawk username (e.g. ``C1234567``)
* Save and close ``run_ppmap.sh"

=================================
Before running PPMAP on a field
=================================

Be sure to replace all instances of ``[USERNAME]`` in the ``run_[FIELD]_all`` file with your Hawk username prior to running PPMAP.
