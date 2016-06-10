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
# Private headers are in src and for CMake need generated config file
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(src)

# Configure as per expat's CMake scripts, always taking defaults
set(XML_CONTEXT_BYTES 1024)
set(XML_DTD 1)
set(XML_NS 1)

include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(TestBigEndian)

check_include_file("dlfcn.h" HAVE_DLFCN_H)
check_include_file("fcntl.h" HAVE_FCNTL_H)
check_include_file("inttypes.h" HAVE_INTTYPES_H)
check_include_file("memory.h" HAVE_MEMORY_H)
check_include_file("stdint.h" HAVE_STDINT_H)
check_include_file("stdlib.h" HAVE_STDLIB_H)
check_include_file("strings.h" HAVE_STRINGS_H)
check_include_file("string.h" HAVE_STRING_H)
check_include_file("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
check_include_file("unistd.h" HAVE_UNISTD_H)

check_function_exists("getpagesize" HAVE_GETPAGESIZE)
check_function_exists("bcopy" HAVE_BCOPY)
check_symbol_exists("memmove" "string.h" HAVE_MEMMOVE)
check_function_exists("mmap" HAVE_MMAP)

#/* Define to 1 if you have the ANSI C header files. */
check_include_files("stdlib.h;stdarg.h;string.h;float.h" STDC_HEADERS)

test_big_endian(WORDS_BIGENDIAN)
#/* 1234 = LIL_ENDIAN, 4321 = BIGENDIAN */
if(WORDS_BIGENDIAN)
    set(BYTEORDER 4321)
else(WORDS_BIGENDIAN)
    set(BYTEORDER 1234)
endif(WORDS_BIGENDIAN)

if(HAVE_SYS_TYPES_H)
    check_symbol_exists("off_t" "sys/types.h" OFF_T)
    check_symbol_exists("size_t" "sys/types.h" SIZE_T)
else(HAVE_SYS_TYPES_H)
    set(OFF_T "long")
    set(SIZE_T "unsigned")
endif(HAVE_SYS_TYPES_H)

configure_file(expat_config.h.cmake expat_config.h)
add_definitions(-DHAVE_EXPAT_CONFIG_H)

if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -wd4996)
endif()

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4expat
    HEADERS
        expat_external.h
        expat.h
    SOURCES
        ${CMAKE_CURRENT_BINARY_DIR}/expat_config.h
        amigaconfig.h
        ascii.h
        asciitab.h
        iasciitab.h
        internal.h
        latin1tab.h
        macconfig.h
        nametab.h
        utf8tab.h
        winconfig.h
        xmlparse.cc
        xmlrole.cc
        xmlrole.h
        xmltok.cc
        xmltok.h
        xmltok_impl.cc
        xmltok_impl.h
        xmltok_ns.cc
)

# List any source specific properties here

