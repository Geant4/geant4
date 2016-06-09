#---Add all subdirectories that contains a CMakeLists.txt file -------------------------------------
file(GLOB_RECURSE files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CMakeLists.txt)
foreach( file ${files} )
  get_filename_component(path ${file} PATH)
  if(path)
    add_subdirectory(${path})
  endif()
endforeach()
