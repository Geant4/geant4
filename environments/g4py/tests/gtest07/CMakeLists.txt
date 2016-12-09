# add teting

add_test(gtest07 python test.py)
configure_file(test.py test.py)
set_property(TEST gtest07
             PROPERTY ENVIRONMENT
             PYTHONPATH=./module:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
           )
