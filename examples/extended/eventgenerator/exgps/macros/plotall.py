#!/usr/bin/env python

from plotfiles import *

for i in range (1,35):
    if i == 23 or i == 27:
        ## skip analysis of problematic macros
        continue
    filename = 'test'
    if i < 10:
        filename = 'test0';
    plot_1_file(filename + str(i))

plot_1_file('test38')

for i in [35, 37]:
    ## skip test36 output producing an error
    plot_2_files('test' + str(i))
