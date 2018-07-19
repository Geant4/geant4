# - Find ZMQ library
#
# It defines:
# ZMQ_FOUND               If the ZeroMQ is found
# ZMQ_INCLUDE_DIRS        PATH to the include directory
# ZMQ_LIBRARIES           libraries
#
# ZMQ_DIR is used for guessing the installed path.
#

message("-- Detecting ZeroMQ package")

find_package(PkgConfig QUIET)
pkg_check_modules(ZMQ QUIET zmq)

if(PC_ZMQ_FOUND)
  set(ZMQ_FOUND 1)
else()
  set(ZMQ_FOUND 0)
endif()

# find include path
find_path(ZMQ_INCLUDE_DIRS zmq.hpp
  HINTS /usr/include ${ZMQ_DIR}/include ${PC_ZMQ_INCLUDE_DIRS})
if(ZMQ_INCLUDE_DIRS MATCHES "^.*-NOTFOUND")
  set(ZMQ_FOUND 0)
endif()

# find library
find_library(ZMQ_LIBRARIES NAMES zmq
  PATHS /usr/lib64 /usr/lib
  ${ZMQ_DIR}/lib64 ${ZMQ_DIR}/lib ${PC_ZMQ_LIBRARY_DIRS})

if(ZMQ_LIBRARIES MATCHES "^.*-NOTFOUND")
  set(ZMQ_FOUND 0)
endif()
#
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(ZMQ DEFAULT_MSG
  ZMQ_INCLUDE_DIRS ZMQ_LIBRARIES)

mark_as_advanced(ZMQ_INCLUDE_DIRS ZMQ_LIBRARY_DIRS ZMQ_LIBRARIES)

message("-- Detecting ZMQ package - done")
