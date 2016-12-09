# add teting

add_test(test12 python test.py)
configure_file(test.py test.py)
set_property(TEST test12 PROPERTY ENVIRONMENT PYTHONPATH=./module)
