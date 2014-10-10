cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
include(CMakeDependentOption)

set(name ChargeExchangeMC)
project(${name})


find_package(Geant4 REQUIRED)

# the project requires Geant4 GDML support!
find_package(gdml QUIET)
if(NOT gdml_FOUND)
  message(STATUS "G4 Examples: GDML package not found. --> ChargeExchangeMC example disabled"
)  
  return()
endif()

# if CEXMC_USE_TCSH is 'yes' then tcsh session support will be built
option(CEXMC_USE_TCSH
    "Build ${name} with tcsh session support
    (requires Geant4 tcsh support)" ON)
if(CEXMC_USE_TCSH)
    find_package(Geant4 COMPONENTS ui_tcsh)
    if(NOT Geant4_ui_tcsh_FOUND)
        message(WARNING
            "Could not find Geant4 tcsh binding, skip tcsh session support")
    endif()
endif()

# if CEXMC_USE_GUI is 'yes' then GUI session support will be built
option(CEXMC_USE_GUI
    "Build ${name} with GUI session support
    (requires Geant4 Qt backend)" ON)
if(CEXMC_USE_GUI)
    find_package(Geant4 COMPONENTS qt)
    if(NOT Geant4_qt_FOUND)
        message(WARNING
            "Could not find Geant4 Qt backend, skip GUI session support")
    endif()
endif()

include(${Geant4_USE_FILE})


set(EXTRA_LIBRARIES )

# if CEXMC_USE_PERSISTENCY is 'yes' then run and events data can be read and
# written; requires boost::serialize headers and library
option(CEXMC_USE_PERSISTENCY
    "Build ${name} with data persistency support
    (requires Boost Serialization library)" ON)
if(CEXMC_USE_PERSISTENCY)
    find_package(Boost COMPONENTS serialization)
    if(Boost_SERIALIZATION_FOUND)
        add_definitions(-DCEXMC_USE_PERSISTENCY)
        include_directories(${Boost_INCLUDE_DIRS})
        link_directories(${Boost_LIBRARY_DIRS})
        list(APPEND EXTRA_LIBRARIES ${Boost_LIBRARIES})
        message(STATUS
            "Library ${Boost_LIBRARIES} was added to the linkage list")
    else()
        message(WARNING
            "Could not find Boost Serialization library, "
            "skip data persistency support")
    endif()
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

# if CEXMC_USE_HISTOGRAMING is 'yes' then ROOT histograming framework will be
# compiled. Notice: if ROOT CERN is not installed in your system then the
# histograming module won't compile anyway
option(CEXMC_USE_HISTOGRAMING
    "Build ${name} with histograming support
    (requires CERN ROOT libraries)" ON)
if(CEXMC_USE_HISTOGRAMING)
    find_package(ROOT QUIET)
    if(ROOT_FOUND)
        add_definitions(-DCEXMC_USE_ROOT)
        include_directories(${ROOT_INCLUDE_DIR})
        list(APPEND EXTRA_LIBRARIES ${ROOT_LIBRARIES})
        message(STATUS
            "Libraries ${ROOT_LIBRARIES} were added to the linkage list")
    else()
        message(WARNING
            "Could not find ROOT package, skip histograming support")
    endif()
endif()

# if CEXMC_USE_ROOTQT is 'yes' then it will be possible to see ROOT histograms
# in GUI session mode
cmake_dependent_option(CEXMC_USE_ROOTQT
    "Build ${name} with GUI histograming support
    (requires CERN ROOT Qt binding)" ON
    "CEXMC_USE_HISTOGRAMING" OFF)
if(CEXMC_USE_ROOTQT AND ROOT_FOUND)
    find_library(ROOTQT_LIBRARY GQt
        PATHS ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARY_DIR}/root
        NO_DEFAULT_PATH)
    if(ROOTQT_LIBRARY)
        add_definitions(-DCEXMC_USE_ROOTQT)
        list(APPEND EXTRA_LIBRARIES -lGQt)
        message(STATUS "Library -lGQt was added to the linkage list")
    else()
        message(WARNING
            "Could not find ROOT Qt library, skip GUI histograming support")
    endif()
endif()

# if CEXMC_USE_GENBOD is 'yes' then original FORTRAN routine GENBOD() will be
# used as phase space generator
option(CEXMC_USE_GENBOD
    "Link ${name} with original FORTRAN function GENBOD()
    (by default C++ reimplementation of GENBOD() is used)" OFF)
if(CEXMC_USE_GENBOD)
    find_program(CERNLIB_CONFIG_EXE cernlib
        PATHS $ENV{CERN_ROOT}/bin)
    if(CERNLIB_CONFIG_EXE)
        execute_process(
            COMMAND ${CERNLIB_CONFIG_EXE} geant321 phtools packlib kernlib
            OUTPUT_VARIABLE CERNLIB_LIBRARIES 
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        add_definitions(-DCEXMC_USE_GENBOD)
        list(APPEND EXTRA_LIBRARIES ${CERNLIB_LIBRARIES})
        message(STATUS
            "Libraries ${CERNLIB_LIBRARIES} were added to the linkage list")
    else()
        message(WARNING
            "Could not find cernlib, skip original FORTRAN GENBOD() support")
    endif()
endif()

# all classes in this project have same member names as their formal constructor
# parameters, so -Wno-shadow option is mandatory to prevent multiple warnings!
add_definitions(-Wno-shadow)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(${name} ${name}.cc ${sources})
target_link_libraries(${name} ${Geant4_LIBRARIES} ${EXTRA_LIBRARIES})

install(TARGETS ${name} DESTINATION bin)

