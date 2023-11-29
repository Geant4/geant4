# Targets file for Build tree. Simply forwards to per-target files if required.
if(PTL_shared_FOUND AND NOT TARGET PTL::ptl-shared)
    include("${CMAKE_CURRENT_LIST_DIR}/ptl-shared.cmake")
endif()

if(PTL_static_FOUND AND NOT TARGET PTL::ptl-static)
    include("${CMAKE_CURRENT_LIST_DIR}/ptl-static.cmake")
endif()
