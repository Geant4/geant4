#------------------------------------------------------------------------------
# sources.cmake
# Module : G4expat
# Package: Geant4.src.G4externals.G4expat
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Created on : 02/06/2011
#
#
#------------------------------------------------------------------------------

# List external includes needed.

# List internal includes needed.
# Private headers are in src!
include_directories(src)

# Add needed compile definitions
if(UNIX AND NOT MSVC)
    add_definitions(-DHAVE_EXPAT_CONFIG_H)
endif()

if(MSVC)
    add_definitions(
        -DWIN32
        -D_WINDOWS
        -D_LIB
        -DCOMPILED_FROM_DSP
        -DXML_UNICODE_WCHAR_T
        -D_VC80_UPGRADE=0x0600
        -D_MBCS)
endif()


#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4expat
    HEADERS
        expat_external.h
        expat.h
        xmlrole.h
        xmltok.h
        xmltok_impl.h
    SOURCES
        amigaconfig.h
        ascii.h
        asciitab.h
        expat_config.h
        expat.h
        iasciitab.h
        internal.h
        latin1tab.h
        macconfig.h
        nametab.h
        utf8tab.h
        winconfig.h
        xmlparse.cc
        xmlrole.cc
        xmltok.cc
        xmltok_impl.c
        xmltok_impl.h
        xmltok_ns.c
)

# List any source specific properties here

