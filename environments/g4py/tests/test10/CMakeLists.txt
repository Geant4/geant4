# add teting

add_test(test10 python test.py)
configure_file(test.py test.py)
set_property(TEST test10 PROPERTY ENVIRONMENT PYTHONPATH=./module)
