# add teting

add_test(test06 python test.py)
configure_file(test.py test.py)
set_property(TEST test06 PROPERTY ENVIRONMENT PYTHONPATH=./module)
