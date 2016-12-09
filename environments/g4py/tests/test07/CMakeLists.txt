# add teting

add_test(test07 python test.py)
configure_file(test.py test.py)
set_property(TEST test07 PROPERTY ENVIRONMENT PYTHONPATH=./module)
