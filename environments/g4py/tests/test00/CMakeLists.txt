# add teting

add_test(test00 python test.py)
configure_file(test.py test.py)
set_property(TEST test00 PROPERTY ENVIRONMENT PYTHONPATH=./module)
