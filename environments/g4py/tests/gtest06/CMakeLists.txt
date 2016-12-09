# add teting

add_test(gtest06 python test.py)
configure_file(test.py test.py)
set_property(TEST gtest06
             PROPERTY ENVIRONMENT
             PYTHONPATH=./module:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
           )
