#
#   Find packages
#

include(FindPackageHandleStandardArgs)
include(MacroUtilities)

ptl_add_interface_library(ptl-external-packages)
ptl_add_interface_library(ptl-threads)
ptl_add_interface_library(ptl-tbb)

################################################################################
#
#                               Threads
#
################################################################################

if(NOT WIN32)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
endif()

find_package(Threads)
if(Threads_FOUND)
    target_link_libraries(ptl-threads INTERFACE Threads::Threads)
    target_link_libraries(ptl-external-packages INTERFACE ptl-threads)
endif()


################################################################################
#
#        TBB
#
################################################################################

if(PTL_USE_TBB)
    find_package(TBB)

    if(TBB_FOUND)
        target_compile_definitions(ptl-tbb INTERFACE PTL_USE_TBB)
        target_include_directories(ptl-tbb SYSTEM INTERFACE ${TBB_INCLUDE_DIRS})
        target_link_libraries(ptl-tbb INTERFACE ${TBB_LIBRARIES})
        target_link_libraries(ptl-external-packages INTERFACE ptl-tbb)
    endif()

endif()
