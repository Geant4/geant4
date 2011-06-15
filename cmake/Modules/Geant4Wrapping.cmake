# - Setup build of any language wrappers around Geant4 C++ libraries
#
# At present the only option here is Python!
# This is yet to be implemented although the general structure is in place
#

#------------------------------------------------------------------------------
# Optional build of Python wrappers (Geant4Py)
# -> Requires Global/Shared library build of Geant4
# -> Requires Boost.Python and Python libraries
#
#include(CMakeDependentOption)
#CMAKE_DEPENDENT_OPTION(GEANT4_WRAP_PYTHON "Build Geant4 Python Wrapping Interface" OFF "NOT GEANT4_BUILD_GRANULAR_LIBS;BUILD_SHARED_LIBS" OFF)
#if(GEANT4_WRAP_PYTHON)
    # We need Boost.Python and Python
    # NB: Watch for FindPythonLibs picking up the static version...
    # Known issue with CMake 2.6 module - see Bug #8389 and 2257 for details
    # of fix
    #find_package(Boost REQUIRED COMPONENTS python)
    #find_package(PythonLibs REQUIRED)
#endif(GEANT4_WRAP_PYTHON)

#GEANT4_ADD_FEATURE(GEANT4_WRAP_PYTHON "Build Geant4 Python Wrapping Interface")


