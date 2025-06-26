#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)

# must be set before project(...) call; version module is needed before
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(Par04)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
SET(GEANT4_MIN_VERSION "11.0")
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 ${GEANT4_MIN_VERSION} REQUIRED ui_all vis_all)
else()
  find_package(Geant4 ${GEANT4_MIN_VERSION} REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup ML-infernce (LCG)
#
option(INFERENCE_LIB "The inference library: currently implemented LWTNN and ONNX" ON)

# LWTNN
if(INFERENCE_LIB)
   find_package(lwtnn QUIET)
   if(lwtnn_FOUND)
	   message("LWTNN inference library found.")
	   add_definitions(-DUSE_INFERENCE)
	   add_definitions(-DUSE_INFERENCE_LWTNN)
   else()
     message("LWTNN not found!")
   endif()
endif()

# ONNX
if(INFERENCE_LIB)
  find_package(OnnxRuntime QUIET)
  find_package(CUDAToolkit QUIET)
   if(OnnxRuntime_FOUND)
       message("ONNX Runtime inference library found.")
   	   add_definitions(-DUSE_INFERENCE)
   	   add_definitions(-DUSE_INFERENCE_ONNX)
      if(CUDAToolkit_FOUND)
        message("Cuda found.")
        add_definitions(-DUSE_CUDA)
      else()
        message("Cuda not found!")
      endif()
   else()
      message("ONNX Runtime not found!")
   endif()
endif()

# TORCH
if(INFERENCE_LIB)
  find_package(Torch QUIET)
   if(Torch_FOUND)
       message("Torch inference library found.")
          add_definitions(-DUSE_INFERENCE)
             add_definitions(-DUSE_INFERENCE_TORCH)
   else()
      message("Torch not found!")
   endif()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(examplePar04 examplePar04.cc ${sources} ${headers})
target_link_libraries(examplePar04 ${Geant4_LIBRARIES} )

if(lwtnn_FOUND)
  target_include_directories(examplePar04 PUBLIC ${lwtnn_INCLUDE_DIRS})
  target_link_libraries(examplePar04 ${lwtnn_LIBRARIES} )
  # Depend on data for runtime
  add_dependencies(examplePar04 examplePar04lwtnndata)
endif()

if(OnnxRuntime_FOUND)
   target_include_directories(examplePar04 PUBLIC ${OnnxRuntime_INCLUDE_DIR})
   target_link_libraries(examplePar04 ${OnnxRuntime_LIBRARY})
    # Cuda_FOUND
    if(CUDAToolkit_FOUND)
      target_link_libraries(examplePar04 CUDA::cudart)     
    endif()
   # Depend on data for runtime
   add_dependencies(examplePar04 examplePar04onnxVAEdata)
   add_dependencies(examplePar04 examplePar04onnxCaloDiTdata)
endif()

if(Torch_FOUND)
   target_include_directories(examplePar04 PUBLIC ${TORCH_INCLUDE_DIRS})
   target_link_libraries(examplePar04 ${TORCH_LIBRARIES})
   message(STATUS "${TORCH_LIBRARIES}")
   # Depend on data for runtime
   add_dependencies(examplePar04 examplePar04torchVAEdata)
   add_dependencies(examplePar04 examplePar04torchCaloDiTdata)
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Par04. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Par04_SCRIPTS
    examplePar04.mac vis.mac common_settings_lowgran.mac common_settings_highgran.mac common_settings_vis.mac common_settings_postInit.mac
  )
if(lwtnn_FOUND)
  set(Par04_SCRIPTS ${Par04_SCRIPTS} examplePar04_lwtnn_vae.mac vis_lwtnn_vae.mac)
endif()
if(OnnxRuntime_FOUND)
  set(Par04_SCRIPTS ${Par04_SCRIPTS} examplePar04_onnx_calodit.mac examplePar04_onnx_vae.mac vis_onnx_calodit.mac vis_onnx_vae.mac)
endif()
if(Torch_FOUND)
  set(Par04_SCRIPTS ${Par04_SCRIPTS} examplePar04_torch_calodit.mac examplePar04_torch_vae.mac vis_torch_calodit.mac vis_torch_vae.mac)
endif()

foreach(_script ${Par04_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Use external data. They are stored outside of the code because of their size.
# THey can be downloaded on demand when this example is built (with inference enabled)
#
include(ExternalProject)
if(lwtnn_FOUND)
  ExternalProject_Add(examplePar04lwtnndata
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/MLModels
    URL https://cern.ch/geant4-data/datasets/examples/extended/parameterisations/Par04/Generator.json
    URL_MD5 af94b1130347428e2e6030930c1ac2cc
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    DOWNLOAD_NO_EXTRACT true
    )
endif()
if(OnnxRuntime_FOUND)
  ExternalProject_Add(examplePar04onnxVAEdata
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/MLModels
    URL https://cern.ch/geant4-data/datasets/examples/extended/parameterisations/Par04/Generator.onnx
    URL_MD5 cacd07c24b704decca28de990850287e
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    DOWNLOAD_NO_EXTRACT true
    )
  ExternalProject_Add(examplePar04onnxCaloDiTdata
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/MLModels
    URL https://cern.ch/geant4-data/datasets/examples/extended/parameterisations/Par04/cd.onnx
    URL_MD5 eb0fa86fc53d9baf72414410a4a9e3c9
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    DOWNLOAD_NO_EXTRACT true
    )
endif()
if(Torch_FOUND)
  ExternalProject_Add(examplePar04torchVAEdata
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/MLModels
    URL https://cern.ch/geant4-data/datasets/examples/extended/parameterisations/Par04/Generator.pt
    URL_MD5 a43337f7f976e976f1127015f2ba61db
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    DOWNLOAD_NO_EXTRACT true
    )
  ExternalProject_Add(examplePar04torchCaloDiTdata
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/MLModels
    URL https://cern.ch/geant4-data/datasets/examples/extended/parameterisations/Par04/cd_cpu.pt
    URL_MD5 c812651390bfc2a4f4a88b3c7f945c56
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    DOWNLOAD_NO_EXTRACT true
    )
endif()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(Par04 DEPENDS examplePar04)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS examplePar04 DESTINATION bin)

