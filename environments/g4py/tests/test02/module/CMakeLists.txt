# - build test library

include_directories(
  ${PYTHON_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
)

link_directories(${Boost_LIBRARY_DIRS})

# library
set(_TARGET test02)
add_library(
  ${_TARGET} SHARED
  AClass.cc
  BClass.cc
  XBase.cc
  test02.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "test02")
set_target_properties(${_TARGET} PROPERTIES SUFFIX ".so")


target_link_libraries (${_TARGET} boost_python ${PYTHON_LIBRARIES})

# install
install(TARGETS ${_TARGET} LIBRARY DESTINATION ${TEST_MODULES_INSTALL_DIR})

