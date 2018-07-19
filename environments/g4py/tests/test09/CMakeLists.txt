# add teting

add_test(test09 python test.py)
configure_file(test.py test.py)
set_property(TEST test09 PROPERTY ENVIRONMENT PYTHONPATH=./module)
