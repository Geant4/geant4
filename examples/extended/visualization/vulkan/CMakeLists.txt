#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(vulkan)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

find_package(Vulkan REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and Vulkan libraries
#
add_executable(exampleVulkan exampleVulkan.cc ${sources} ${headers})
target_link_libraries(exampleVulkan ${Geant4_LIBRARIES} Vulkan::Vulkan)

#----------------------------------------------------------------------------
# Copy all scripts and benchmark files to the build directory
#
set(EXAMPLEVULKAN_SCRIPTS
  init_vis.mac
  vis.mac
  vulkan.mac
  opengl.mac
  benchmark_opengl_vs_vulkan.json
  )

foreach(_script ${EXAMPLEVULKAN_SCRIPTS})
  if(EXISTS ${PROJECT_SOURCE_DIR}/${_script})
    configure_file(
      ${PROJECT_SOURCE_DIR}/${_script}
      ${PROJECT_BINARY_DIR}/${_script}
      COPYONLY
      )
  endif()
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use
#
add_custom_target(vulkan DEPENDS exampleVulkan)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleVulkan DESTINATION bin)

