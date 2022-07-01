# - G4globman module build definition

# - Configure header for preprocessor symbols
# Convert CMake variables -> #cmakedefine symbols
set(G4MULTITHREADED ${GEANT4_BUILD_MULTITHREADED})
set(G4_STORE_TRAJECTORY ${GEANT4_BUILD_STORE_TRAJECTORY})
set(G4VERBOSE ${GEANT4_BUILD_VERBOSE_CODE})

configure_file(${CMAKE_CURRENT_LIST_DIR}/include/G4GlobalConfig.hh.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/G4GlobalConfig.hh)

geant4_get_datasetnames(GEANT4_DATASETS)
foreach(_ds ${GEANT4_DATASETS})
  geant4_get_dataset_property(${_ds} ENVVAR ${_ds}_ENVVAR)
  geant4_get_dataset_property(${_ds} DIRECTORY ${_ds}_DIR)
  set(GEANT4_DATASET_LIST
    "${GEANT4_DATASET_LIST}\n{ \"${${_ds}_ENVVAR}\", \"${${_ds}_DIR}\" },")
endforeach()

configure_file(${CMAKE_CURRENT_LIST_DIR}/include/G4FindDataDir.hh.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/G4FindDataDir.hh)

#
# Define the Geant4 Module.
#
geant4_add_module(G4globman
  PUBLIC_HEADERS
    ${CMAKE_CURRENT_BINARY_DIR}/include/G4GlobalConfig.hh
    G4Allocator.hh
    G4AllocatorList.hh
    G4AllocatorPool.hh
    G4ApplicationState.hh
    G4AutoDelete.hh
    G4AutoLock.hh
    G4Backtrace.hh
    G4BuffercoutDestination.hh
    G4CacheDetails.hh
    G4Cache.hh
    G4coutDestination.hh
    G4coutFormatters.hh
    G4DataVector.hh
    G4DataVector.icc
    G4EnvironmentUtils.hh
    G4ErrorPropagatorData.hh
    G4ErrorPropagatorData.icc
    G4Evaluator.hh
    G4Exception.hh
    G4ExceptionSeverity.hh
    G4Exp.hh
    G4FastVector.hh
    G4FilecoutDestination.hh
    G4Filesystem.hh
    G4FPEDetection.hh
    G4GeometryTolerance.hh
    G4GlobalConfig.hh.in
    G4ios.hh
    G4LockcoutDestination.hh
    G4Log.hh
    G4MasterForwardcoutDestination.hh
    G4MTBarrier.hh
    G4MTcoutDestination.hh
    G4MulticoutDestination.hh
    G4OrderedTable.hh
    G4PhysicalConstants.hh
    G4Physics2DVector.hh
    G4Physics2DVector.icc
    G4PhysicsFreeVector.hh
    G4PhysicsLinearVector.hh
    G4PhysicsLogVector.hh
    G4PhysicsModelCatalog.hh
    G4PhysicsOrderedFreeVector.hh
    G4PhysicsTable.hh
    G4PhysicsTable.icc
    G4PhysicsVector.hh
    G4PhysicsVector.icc
    G4PhysicsVectorType.hh
    G4Pow.hh
    G4Profiler.hh
    G4Profiler.icc
    G4ReferenceCountedHandle.hh
    G4RotationMatrix.hh
    G4SIunits.hh
    G4SliceTimer.hh
    G4SliceTimer.icc
    G4StateManager.hh
    G4StateManager.icc
    G4String.hh
    G4String.icc
    G4strstreambuf.hh
    G4strstreambuf.icc
    G4SystemOfUnits.hh
    G4TaskGroup.hh
    G4Task.hh
    G4TaskManager.hh
    G4TaskSingletonDelegator.hh
    G4TBBTaskGroup.hh
    G4ThreadData.hh
    G4Threading.hh
    G4ThreadLocalSingleton.hh
    G4ThreadPool.hh
    G4ThreeVector.hh
    G4TiMemory.hh
    G4Timer.hh
    G4Timer.icc
    G4Tokenizer.hh
    G4TWorkspacePool.hh
    G4TwoVector.hh
    G4Types.hh
    G4UnitsTable.hh
    G4UnitsTable.icc
    G4UserLimits.hh
    G4UserLimits.icc
    G4UserTaskQueue.hh
    G4Version.hh
    G4VExceptionHandler.hh
    G4VNotifier.hh
    G4VStateDependent.hh
    G4VTask.hh
    G4VUserTaskQueue.hh
    globals.hh
    templates.hh
    tls.hh
    windefs.hh
  SOURCES
    G4Allocator.cc
    G4AllocatorPool.cc
    G4AllocatorList.cc
    G4BuffercoutDestination.cc
    G4CacheDetails.cc
    G4coutDestination.cc
    G4coutFormatters.cc
    G4DataVector.cc
    G4ErrorPropagatorData.cc
    G4Exception.cc
    G4FilecoutDestination.cc
    G4FindDataDir.cc
    G4GeometryTolerance.cc
    G4ios.cc
    G4LockcoutDestination.cc
    G4MasterForwardcoutDestination.cc
    G4MTBarrier.cc
    G4MTcoutDestination.cc
    G4OrderedTable.cc
    G4PhysicsFreeVector.cc
    G4PhysicsLinearVector.cc
    G4PhysicsLogVector.cc
    G4PhysicsModelCatalog.cc
    G4PhysicsTable.cc
    G4PhysicsVector.cc
    G4Physics2DVector.cc
    G4Pow.cc
    G4Profiler.cc
    G4ReferenceCountedHandle.cc
    G4SliceTimer.cc
    G4StateManager.cc
    G4ThreadLocalSingleton.cc
    G4Threading.cc
    G4Timer.cc
    G4UnitsTable.cc
    G4VExceptionHandler.cc
    G4VStateDependent.cc)

# - Add path to generated header
geant4_module_include_directories(G4globman
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)

geant4_module_link_libraries(G4globman
  PUBLIC
    ${CLHEP_LIBRARIES}
    ${timemory_LIBRARIES}
    ${PTL_LIBRARIES}
    ${GEANT4_CXX_FILESYSTEM_LIBRARY} # to temporarily support libstdc++fs, libc++fs
    )

# TEMP WORKAROUND: When building/testing examples uing ROOT, ROOT's
# dictionary generation is not smart enough to handle target usage
# requirements for include paths. Explicitly add the path to the
# generated header into build time include paths...
# We eventually want to do this through "..._include_directories" or to
# remove entirely and require a min version of ROOT
set_property(GLOBAL APPEND
  PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/include")
