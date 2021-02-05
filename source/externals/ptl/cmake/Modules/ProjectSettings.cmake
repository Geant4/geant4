################################################################################
#
# Project settings
#
################################################################################

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)

add_feature(CMAKE_C_FLAGS_${_CONFIG} "C compiler build-specific flags")
add_feature(CMAKE_CXX_FLAGS_${_CONFIG} "C++ compiler build-specific flags")

################################################################################
#
#   installation directories
#
################################################################################

# cmake installation folder
if(NOT CMAKE_INSTALL_CONFIGDIR)
  set(CMAKE_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})
endif()


# create the full path version and generic path versions
foreach(_TYPE in DATAROOT CMAKE INCLUDE LIB BIN MAN DOC)
    # set the absolute versions
    if(NOT IS_ABSOLUTE "${CMAKE_INSTALL_${_TYPE}DIR}")
        set(CMAKE_INSTALL_FULL_${_TYPE}DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_${_TYPE}DIR})
    else(NOT IS_ABSOLUTE "${CMAKE_INSTALL_${_TYPE}DIR}")
        set(CMAKE_INSTALL_FULL_${_TYPE}DIR ${CMAKE_INSTALL_${_TYPE}DIR})
    endif(NOT IS_ABSOLUTE "${CMAKE_INSTALL_${_TYPE}DIR}")

    # generic "PROJECT_INSTALL_" variables (used by documentation)"
    set(PROJECT_INSTALL_${_TYPE}DIR ${CMAKE_INSTALL_${TYPE}DIR})
    set(PROJECT_INSTALL_FULL_${_TYPE}DIR ${CMAKE_INSTALL_FULL_${TYPE}DIR})

endforeach(_TYPE in DATAROOT CMAKE INCLUDE LIB BIN MAN DOC)
