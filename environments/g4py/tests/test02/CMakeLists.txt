# add teting

add_test(test02 python test.py)
configure_file(test.py test.py)
set_property(TEST test02 PROPERTY ENVIRONMENT PYTHONPATH=./module)
