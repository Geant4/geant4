# - Build for Geant4-vendored version of Expat
#
# It largely reproduces Expat's core CMake script for checks/defines
# but replaces the library build with Geant4 specifics for naming
# an installing.

#----------------------------------------------------------------------------
# - Expat core build
set(XML_CONTEXT_BYTES 1024)
set(XML_DTD 1)
set(XML_NS 1)

include(CheckCSourceCompiles)
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

check_symbol_exists("getpagesize" "unistd.h" HAVE_GETPAGESIZE)
check_symbol_exists("mmap" "sys/mman.h" HAVE_MMAP)
check_symbol_exists("getrandom" "sys/random.h" HAVE_GETRANDOM)

check_symbol_exists("arc4random_buf" "${_bsd}stdlib.h" HAVE_ARC4RANDOM_BUF)
if(NOT HAVE_ARC4RANDOM_BUF)
    check_symbol_exists("arc4random" "${_bsd}stdlib.h" HAVE_ARC4RANDOM)
endif()

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

check_c_source_compiles("
        #include <stdlib.h>  /* for NULL */
        #include <unistd.h>  /* for syscall */
        #include <sys/syscall.h>  /* for SYS_getrandom */
        int main() {
            syscall(SYS_getrandom, NULL, 0, 0);
            return 0;
        }"
    HAVE_SYSCALL_GETRANDOM)

configure_file(expat_config.h.cmake expat_config.h)
add_definitions(-DHAVE_EXPAT_CONFIG_H)

if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -wd4996)
endif()

#----------------------------------------------------------------------------
# - Geant4 Integration
# Use geant4_add_library to create G4expat and static variant if needed
set(G4expat_HEADERS
  include/expat_external.h
  include/expat.h
  ${CMAKE_CURRENT_BINARY_DIR}/expat_config.h)

set(G4expat_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/expat_config.h
  src/ascii.h
  src/asciitab.h
  src/iasciitab.h
  src/internal.h
  src/latin1tab.h
  src/nametab.h
  src/siphash.h
  src/utf8tab.h
  src/winconfig.h
  src/xmlparse.c
  src/xmlrole.c
  src/xmlrole.h
  src/xmltok.c
  src/xmltok.h
  src/xmltok_impl.c
  src/xmltok_impl.h
  src/xmltok_ns.c)

# - Create target(s) and set usage requirements
geant4_add_external_library(NAME G4expat SOURCES ${G4expat_SOURCES} ${G4expat_HEADERS})
foreach(__g4expat_target G4expat G4expat-static)
  if(TARGET ${__g4expat_target})
    target_include_directories(${__g4expat_target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    )
  endif()
endforeach()

if(TARGET G4expat-static)
  target_compile_definitions(G4expat-static PUBLIC XML_STATIC)
endif()

# Install headers alongside Geant4
install(FILES ${G4expat_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  COMPONENT Development)
