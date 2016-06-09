# - Geant4MakeRules_cxx
# Sets the default make rules for a CXX build, specifically the initialization
# of the compiler flags on a platform dependent basis
#

message(STATUS "setting default compiler flags for CXX")


#------------------------------------------------------------------------------
# GNU C++ Compiler on all(?) platforms
#
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_INIT "-W -Wall -ansi -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -pipe")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g") 
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
endif()


#------------------------------------------------------------------------------
# MSVC - all (?) versions
# I don't understand VS flags at all.... Guess try the CMake defaults first
# and see what happens!
if(MSVC)
    # Hmm, WIN32-VC.gmk uses dashes, but cmake uses slashes, latter probably
    # best for native build.
    set(CMAKE_CXX_FLAGS_INIT "-MD -GR -EHsc -Zm200 -nologo -D_CONSOLE -D_WIN32 -DWIN32 -DOS")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Od -Zi") 
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

    # We may also have to set linker flags....
endif()


#------------------------------------------------------------------------------
# Intel C++ Compilers - all (?) platforms
#
# Sufficient id on all platforms?
if(CMAKE_CXX_COMPILER MATCHES "icpc.*")
    set(CMAKE_CXX_FLAGS_INIT "-w1 -Wnon-virtual-dtor -Wpointer-arith -Wwrite-strings -ansi -fp-model precise")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g") 
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
endif()


#------------------------------------------------------------------------------
# Ye Olde *NIX/Compiler Systems
# NB: *NOT* Supported... Only provided as legacy.
# None are tested...
# Whilst these use flags taken from existing Geant4 setup, may want to see if
# CMake defaults on these platforms are good enough...
#
if(UNIX AND NOT CMAKE_COMPILER_IS_GNUCXX)
    #--------------------------------------------------------------------------
    # IBM xlC compiler
    #
    if(CMAKE_CXX_COMPILER MATCHES "xlC")
        set(CMAKE_CXX_FLAGS_INIT "")
        set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g -qdbextra -qcheck=all -qfullpath -qtwolink -+") 
        set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -qtwolink -+")
        set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O3 -qtwolink -+")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O3 -g -qdbextra -qcheck=all -qfullpath -qtwolink -+")
    endif()

    #--------------------------------------------------------------------------
    # HP aC++ Compiler
    #
    if(CMAKE_CXX_COMPILER MATCHES "aCC")
        set(CMAKE_CXX_FLAGS_INIT "+DAportable +W823")
        set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g") 
        set(CMAKE_CXX_FLAGS_RELEASE_INIT "+O2 +Onolimit")
        set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O2 +Onolimit")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 +Onolimit -g")
    endif()

    #--------------------------------------------------------------------------
    # IRIX MIPSpro CC Compiler
    #
    if(CMAKE_CXX_COMPILER MATCHES "CC" AND CMAKE_SYSTEM_NAME MATCHES "IRIX")
        set(CMAKE_CXX_FLAGS_INIT "-ptused -DSOCKET_IRIX_SOLARIS")
        set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g") 
        set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O -OPT:Olimit=5000")
        set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O -OPT:Olimit=5000")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O -OPT:Olimit=5000 -g")
    endif()

    #--------------------------------------------------------------------------
    # SunOS CC Compiler
    # - CMake may do a reasonable job on its own here...
endif()

