# - add tests components

set(TEST_MODULES_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/tests)

# tests for boost_python
add_subdirectory(test00/module)
add_subdirectory(test01/module)
add_subdirectory(test02/module)
#add_subdirectory(test03/module EXCLUDE_FROM_ALL)
add_subdirectory(test04/module)
add_subdirectory(test05/module)
add_subdirectory(test06/module)
add_subdirectory(test07/module)
add_subdirectory(test08/module)
add_subdirectory(test09/module)
add_subdirectory(test10/module)
add_subdirectory(test11/module)
add_subdirectory(test12/module)
add_subdirectory(test13/module)

# tests for g4py
add_subdirectory(gtest01/module)
