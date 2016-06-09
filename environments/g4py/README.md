Geant4Py
========

_A set of python modules for using Geant4_


System Requirements
-------------------
### Python

Python2.x and Python 3.x (experimental)

### Geant4

built with CMake

### CLHEP

use Geant4-CLHEP

### (Optional)

ROOT for histogramming/analysis


How to Install
--------------
There is a configuration script for building the package.

    # ./configure --help
    `configure' configures Geant4Py to adapt to many kinds of systems.
    
    Usage:  ./configure SYSTEM [OPTION]... [VAR=VALUE]...

    SYSTEM: System type (see Supported Arhitectures)

    Options:
      -h, --help                Display this help and exit
    
    Installation directories:
      --prefix=PREFIX           Installation prefix  [./]
      --libdir=DIR              Python modules dir [PREFIX/lib]

    Fine tuning of the library path:
      --with-g4install-dir=DIR  Geant4 installed dir
 
      --with-python-incdir=DIR  Python header dir [/usr/include/python(2.#)],
                                (location of pyconfig.h)
      --with-python-libdir=DIR  Python library dir [/usr/lib(64)]
      --with-python3            Use Python3
    
     --with-boost-incdir=DIR   BOOST-C++ header dir [/usr/include],
                                (location of boost/)
     --with-boost-libdir=DIR   BOOST-C++ library dir [/usr/lib]
     --with-boost-python-lib=LIB library name of libboost_python.so [boost_python]

     --with-extra-dir=DIR      Install path for extra packages [/usr/local]

     --with-xercesc-incdir=DIR Xerces-C header dir [/usr/include]
     --with-xercesc-libdir=DIR Xerces-C library dir [/usr/lib(64)]

    Enable/disable options: prefix with either --enable- or --disable-
      openglx      OpenGLX support    [auto]
      openglxm     OpenGLXm support   [disable]
      raytracerx   RayTracerX support [disable]

    Supported Architectures:
      linux           for Linux gcc (32bit)
      linux64         for Linux gcc (64bit)
      linuxx8664gcc   for Linux gcc (64bit)
      macosx          for Apple OS X with gcc


For example,

    # ./configure linux64 -with-g4install-dir=<geant4 install path with CMake>

The configuration script will create config/config.gmk, 
which describes your environment.

After executing the `configure` script, then 

    # make
    # make install


How to Use:
-----------
Some environment variables are required at run time.

### *PYTHONPATH*

Python module search directories, given by a colon-separated list of directories. 
Practically, `$(INSTALL_PATH)/lib` is ok.

### *LD\_LIBRARY\_PATH*

Additional shared library path to be searched.

Some libraries paths are already specified via *-rpath* linker option, 
so these paths do not have to be added to *LD\_LIBRARY\_PATH*.

* Geant4
* boost-python library
* $(INSTALL_PATH)/lib
* $(G4PY\_EXLIB\_LIBDIR) specified in a module makefile.


You can import Geant4Py modules in Python just like

    >>> import Geant4

