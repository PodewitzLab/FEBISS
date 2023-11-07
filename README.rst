Free Energy Based Identification of Solvation Sites (FEBISS)
============================================================

Installation
------------

FEBISS can be installed using pip (pip3) once the repository has been cloned:

.. code-block:: bash

   git clone https://github.com/PodewitzLab/FEBISS.git
   pip install ./FEBISS/

A non super user can install the package using a virtual environment, or
the ``--user`` flag.

A manual with detailed instructions can be found in the Github repo.


Prerequisites
-------------

Basic Requirements
..................

FEBISS is expected to run on Linux systems and the following
programs/packages are required:

- Python3
- Git
- GCC >= v7.0.0

Febiss Python Package
.....................

The main Python package called FEBISS requires only basic additional packages, which will 
be automatically installed when installing FEBISS, using ``pip``. These packages are
listed in the file ``requirements.txt``.

C++ Requirements
................

To run analyses of trajectories the open-source software CPPTRAJ modified with the GIGIST 
repository is used. These dependencies are not necessary for the installation of FEBISS, 
but rather FEBISS provides a script to set-up these dependencies via

.. code-block:: bash

   febiss_setup

Please be aware that CPPTRAJ may require libraries that cannot be installed via FEBISS, but 
have to be installed by the user first. CPPTRAJ makes use of the following libraries:

- NetCDF
- BLAS
- LAPACK
- Gzip
- Bzip2
- Parallel NetCDF (-mpi build only, for NetCDF trajectory output in parallel)
- CUDA (-cuda build only)
- FFTW (mostly optional; required for PME functionality and very large FFTs)

We therefore recommend to install some basic libraries via

.. code-block:: bash

   sudo apt-get install libblas-dev liblapack-dev libbz2-dev libnetcdf-dev

Should you encounter difficulties in the installation of CPPTRAJ, we refer to the README
of the GIGIST and CPPTRAJ repositories.

Basic Usage
-----------

If you installed FEBISS and its dependencies ``febiss_setup``, you can use FEBISS 
to analyse trajectories for water placements, plot the data and select the waters you want 
to further investigate within a bar chart.

To get a list of all available options and a useful input file, you can call

.. code-block:: bash

   febiss_settings

This will place the file ``all-settings.yaml`` in your current directory. This input 
file requires only 2 alterations to be a valid input file. You have to give the name of your 
topology file and the base name of your trajectory file(s). Along those 2 mandatory settings 
you find all other available settings both for the GIST analysis and the plotting of the retrieved 
data. Once you performed the analysis, you can also skip this step and directly plot the data.
The analysis and plotting are done via calling the main program along with the yaml input file:

.. code-block:: bash

   febiss all-settings.yaml


