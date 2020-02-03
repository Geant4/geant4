# add teting

add_test(test00 python test.py)

if (${PYTHON_VERSION_MAJOR} MATCHES "2")
  configure_file(test.py test.py)
elseif (${PYTHON_VERSION_MAJOR} MATCHES "3")
  configure_file(python3/test.py test.py)
endif()

set_property(TEST test00 PROPERTY ENVIRONMENT PYTHONPATH=./module)
