Geant4Py
========

_A set of python modules for using Geant4_


System Requirements
-------------------

### CMake
Building system is migrated to CMake system.

### Python

Python2.x and Python 3.x (experimental)

### Boost

Boost_Python is needed.

### (Optional)

ROOT for histogramming/analysis

### Tested Platforms

* CentOS6 (RH6 clones)
* OSX Yosemite


How to Install
--------------
Before building library, *GEANT4_INSTALL* environment variable should be set.

    # export GEANT4_INSTALL=<Geant4 install path> (zsh, bash)
    # setenv GEANT4_INSTALL <Geant4 install path> (csh)
    (G4 install path is the path specified by "CMAKE_INSTALL_PREFIX" when building Geant4) 
 
then

    # mkdir build
    # cd build
    # cmake ..   
    # make
    # make install

If you want to run the testing component,

    # cd build/tests
    # make; make install

By default, g4py is installed in \<g4py\>/lib(64) directory.


How to Use:
-----------
Some environment variables are required at run time.

### *PYTHONPATH*

Python module search directories, given by a colon-separated list of directories, like


    # export PYTHONPATH=<g4py>/lib64:<g4py>/examples:<g4py>/tests  (zsh, bash)
    # setenv PYTHONPATH <g4py>/lib64:<g4py>/examples:<g4py>/tests (csh)



### Getting started
You can import Geant4Py modules in Python just like

    >>> import Geant4

