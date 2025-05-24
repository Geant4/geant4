add_library(G4cpp-httplib INTERFACE)

target_include_directories(G4cpp-httplib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)