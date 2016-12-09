# add teting

add_test(test11 python test.py)
configure_file(test.py test.py)
set_property(TEST test11 PROPERTY ENVIRONMENT PYTHONPATH=./module)
