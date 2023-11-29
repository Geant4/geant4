# - G4zlib build definition

# -----
# Primary ZLIB CMake script, adapted

include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckCSourceCompiles)

check_include_file(sys/types.h HAVE_SYS_TYPES_H)
check_include_file(stdint.h    HAVE_STDINT_H)
check_include_file(stddef.h    HAVE_STDDEF_H)

#
# Check to see if we have large file support
#
set(CMAKE_REQUIRED_DEFINITIONS -D_LARGEFILE64_SOURCE=1)
# We add these other definitions here because CheckTypeSize.cmake
# in CMake 2.4.x does not automatically do so and we want
# compatibility with CMake 2.4.x.
if(HAVE_SYS_TYPES_H)
    list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_SYS_TYPES_H)
endif()
if(HAVE_STDINT_H)
    list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_STDINT_H)
endif()
if(HAVE_STDDEF_H)
    list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_STDDEF_H)
endif()
check_type_size(off64_t OFF64_T)
if(HAVE_OFF64_T)
   add_definitions(-D_LARGEFILE64_SOURCE=1)
endif()
set(CMAKE_REQUIRED_DEFINITIONS) # clear variable

#
# Check for fseeko
#
check_function_exists(fseeko HAVE_FSEEKO)
if(NOT HAVE_FSEEKO)
    add_definitions(-DNO_FSEEKO)
endif()

#
# Check for unistd.h
#
check_include_file(unistd.h Z_HAVE_UNISTD_H)

if(MSVC)
#    set(CMAKE_DEBUG_POSTFIX "d")
    add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
    add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
	add_definitions(-DZLIB_DLL)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    # If we're doing an out of source build and the user has a zconf.h
    # in their source tree...
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/zconf.h)
        message(STATUS "Renaming")
        message(STATUS "    ${CMAKE_CURRENT_SOURCE_DIR}/zconf.h")
        message(STATUS "to 'zconf.h.included' because this file is included with zlib")
        message(STATUS "but CMake generates it automatically in the build directory.")
        file(RENAME ${CMAKE_CURRENT_SOURCE_DIR}/zconf.h ${CMAKE_CURRENT_SOURCE_DIR}/zconf.h.included)
  endif()
endif()

configure_file(	${CMAKE_CURRENT_SOURCE_DIR}/src/zconf.h.cmakein
		${CMAKE_CURRENT_BINARY_DIR}/zconf.h @ONLY)

# -----
# - Geant4 specific part to integrate
#
# Intel/Clang may warn about -Wdeprecated-non-prototype, but per https://github.com/madler/zlib/issues/633
# we suppress this warning if the compiler supports the flag
include(CheckCCompilerFlag)
check_c_compiler_flag("-Wno-deprecated-non-prototype" G4ZLIB_NEEDS_DNP)

# Headers listed under Sources are internal zlib headers
# Private headers are in src!
set(ZLIB_PUBLIC_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/zconf.h
  include/zlib.h )
set(ZLIB_PRIVATE_HDRS
  include/crc32.h
  include/deflate.h
  include/gzguts.h
  include/inffast.h
  include/inffixed.h
  include/inflate.h
  include/inftrees.h
  include/trees.h
  include/zutil.h )
set(ZLIB_SRCS
  src/adler32.c
  src/compress.c
  src/crc32.c
  src/deflate.c
  src/gzclose.c
  src/gzlib.c
  src/gzread.c
  src/gzwrite.c
  src/inflate.c
  src/infback.c
  src/inftrees.c
  src/inffast.c
  src/trees.c
  src/uncompr.c
  src/zutil.c )

geant4_add_external_library(NAME G4zlib SOURCES ${ZLIB_PUBLIC_HDRS} ${ZLIB_PRIVATE_HDRS} ${ZLIB_SRCS})

# Add header paths for private/generated files in Build, plus flags if present
foreach(__g4zlib_target G4zlib G4zlib-static)
  if(TARGET ${__g4zlib_target})
    target_include_directories(${__g4zlib_target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    )

    if(G4ZLIB_NEEDS_DNP)
      target_compile_options(${__g4zlib_target} PRIVATE "-Wno-deprecated-non-prototype")
    endif()
  endif()
endforeach()

# Install headers alongside Geant4
install(FILES ${ZLIB_PUBLIC_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  COMPONENT Development)
