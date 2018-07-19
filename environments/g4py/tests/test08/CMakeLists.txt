# add teting

add_test(test08 python test.py)
configure_file(test.py test.py)
set_property(TEST test08 PROPERTY ENVIRONMENT PYTHONPATH=./module)
