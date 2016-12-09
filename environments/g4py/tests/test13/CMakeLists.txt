# add teting

add_test(test13 python test.py)
configure_file(test.py test.py)
set_property(TEST test13 PROPERTY ENVIRONMENT PYTHONPATH=./module)
