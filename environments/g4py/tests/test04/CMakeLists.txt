# add teting

add_test(test04-1 python test1.py)
configure_file(test1.py test1.py)
set_property(TEST test04-1 PROPERTY ENVIRONMENT PYTHONPATH=./module)
#
add_test(test04-2 python test2.py)
configure_file(test2.py test2.py)
set_property(TEST test04-2 PROPERTY ENVIRONMENT PYTHONPATH=./module)
