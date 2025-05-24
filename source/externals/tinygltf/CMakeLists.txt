add_library(G4tinygltf INTERFACE)

target_include_directories(G4tinygltf INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)