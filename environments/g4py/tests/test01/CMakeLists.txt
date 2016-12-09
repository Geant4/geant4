# add teting

add_test(test01 python test.py)
configure_file(test.py test.py)
set_property(TEST test01 PROPERTY ENVIRONMENT PYTHONPATH=./module)
