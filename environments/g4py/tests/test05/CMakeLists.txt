# add teting

add_test(test05 python test.py)
configure_file(test.py test.py)
set_property(TEST test05 PROPERTY ENVIRONMENT PYTHONPATH=./module)
