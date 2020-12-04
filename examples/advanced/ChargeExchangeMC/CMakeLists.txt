cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
include(CMakeDependentOption)

set(name ChargeExchangeMC)
project(${name})

find_package(Geant4 REQUIRED gdml)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

# Setup of ROOT ANALYSIS : optional.
find_package(ROOT QUIET)
cmake_dependent_option(CEXMC_USE_HISTOGRAMING "Build example with analysis objects" ON "ROOT_FOUND" OFF)

set(EXTRA_LIBRARIES )

# if CEXMC_USE_PERSISTENCY is 'yes' then run and events data can be read and
# written; requires boost::serialize headers and library
option(CEXMC_USE_PERSISTENCY
    "Build ${name} with data persistency support
    (requires Boost Serialization library)" OFF)

if(CEXMC_USE_PERSISTENCY)
   find_package(Boost REQUIRED serialization)
   add_definitions(-DCEXMC_USE_PERSISTENCY)
   list(APPEND EXTRA_LIBRARIES Boost::serialization)
    message(STATUS "Library Boost::serialization was added to the linkage list")
endif()

# if CEXMC_USE_CUSTOM_FILTER is 'yes' then Custom filter can be used for
# existing events data; requires boost::spirit 2.x headers. Notice: if
# CEXMC_USE_PERSISTENCY is not 'yes' then Custom Filter will not be used anyway
cmake_dependent_option(CEXMC_USE_CUSTOM_FILTER
    "Build ${name} with custom filter support
    (requires Boost Spirit library)" OFF
    "CEXMC_USE_PERSISTENCY" OFF)
if(CEXMC_USE_CUSTOM_FILTER)
    add_definitions(-DCEXMC_USE_CUSTOM_FILTER)
endif()

# if CEXMC_DEBUG_CUSTOM_FILTER is 'yes' then AST trees will be printed out
cmake_dependent_option(CEXMC_DEBUG_CUSTOM_FILTER
    "Debug custom filter" OFF
    "CEXMC_USE_CUSTOM_FILTER" OFF)
if(CEXMC_USE_CUSTOM_FILTER)
    add_definitions(-DCEXMC_DEBUG_CF)
endif()

# if CEXMC_USE_QGSP_BERT is 'yes' then QGSP_BERT will be used as basic physics,
# otherwise - FTFP_BERT or QGSP_BIC_EMY
option(CEXMC_USE_QGSP_BERT
    "Build ${name} with QGSP_BERT physics list
    (default physics list is FTFP_BERT)" OFF)
if(CEXMC_USE_QGSP_BERT)
        add_definitions(-DCEXMC_USE_QGSP_BERT)
endif()

# if CEXMC_USE_QGSP_BIC_EMY is 'yes' then QGSP_BIC_EMY will be used as basic
# physics, otherwise - FTFP_BERT or QGSP_BERT
cmake_dependent_option(CEXMC_USE_QGSP_BIC_EMY
    "Build ${name} with QGSP_BIC_EMY physics list
    (default physics list is FTFP_BERT)" OFF
    "NOT CEXMC_USE_QGSP_BERT" OFF)
if(CEXMC_USE_QGSP_BIC_EMY)
    add_definitions(-DCEXMC_USE_QGSP_BIC_EMY)
endif()

# if CEXMC_DEBUG_TP is 'yes' then additional info will be printed on track
# points data
option(CEXMC_DEBUG_TP
    "Print debug information for track points" OFF)
if(CEXMC_DEBUG_TP)
    add_definitions(-DCEXMC_DEBUG_TP)
endif()

if(CEXMC_USE_HISTOGRAMING)
	EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
	EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
endif(CEXMC_USE_HISTOGRAMING)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(${name} ${name}.cc ${sources})

target_link_libraries(${name} ${Geant4_LIBRARIES} ${EXTRA_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ChargeExchangeMC. This is so that we can run the executable
# directly because it relies on these scripts being in the current working
# directory.

set(chargeexchangemc_SCRIPTS
    batch.mac replay.mac init.mac preinit.mac lht.gdml)

foreach(_script ${chargeexchangemc_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${name} DESTINATION bin)

