
# include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found.")  
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${ROOTSYS}/include)

### file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
### file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME CommonSW
    HEADERS
       Beam.hh
       ExecBase.hh
       ExecProcessLevel.hh
       FTFPWrapper.hh
       ProcessWrapper.hh
       QGSPWrapper.hh
       TstDiscreteProcessReader.hh
       TstHisto.hh
       TstHistoSet.hh
       TstPhysListReader.hh
       TstPrimaryGeneratorAction.hh
       TstReader.hh
       TstTarget.hh
    SOURCES
       Beam.cc
       ExecBase.cc
       ExecProcessLevel.cc
       FTFPWrapper.cc
       ProcessWrapper.cc
       QGSPWrapper.cc
       TstDiscreteProcessReader.cc
       TstHisto.cc
       TstHistoSet.cc
       TstPhysListReader.cc
       TstPrimaryGeneratorAction.cc
       TstReader.cc
       TstTarget.cc
    GLOBAL_DEPENDENCIES
       ${Geant4_LIBRARIES}
       ${ROOT_LIBRARIES}
    LINK_LIBRARIES
)
