# -------------------------------------------------------------------------------------- #
#
# installation directories
#
# -------------------------------------------------------------------------------------- #
include(GNUInstallDirs)

# cmake installation folder (not exposed to user)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
if(NOT IS_ABSOLUTE "${CMAKE_INSTALL_CMAKEDIR}")
    set(CMAKE_INSTALL_FULL_CMAKEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_CMAKEDIR})
else()
    set(CMAKE_INSTALL_FULL_CMAKEDIR ${CMAKE_INSTALL_CMAKEDIR})
endif()

# Create PTL_ versions that are used internally and may be overriden by master projects
# When we are the master project, these are *always* set to the value(s) of
# CMAKE_INSTALL_<TYPE> Otherwise, the value of CMAKE_INSTALL_<TYPE> is used as the
# default.
foreach(_TYPE DATAROOT CMAKE INCLUDE LIB BIN MAN DOC)
    if(PTL_MASTER_PROJECT)
        set(PTL_INSTALL_${_TYPE}DIR "${CMAKE_INSTALL_${_TYPE}DIR}")
    else()
        if(NOT DEFINED PTL_INSTALL_${_TYPE}DIR)
            set(PTL_INSTALL_${_TYPE}DIR "${CMAKE_INSTALL_${_TYPE}DIR}")
        endif()
    endif()
endforeach()
