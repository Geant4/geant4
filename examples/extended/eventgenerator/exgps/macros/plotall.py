#!/usr/bin/python

from plotfiles import *

for i in range (1,35):
    plot_1_file("test"+str(i))

plot_1_file('test38')

for i in range (35,38):
    plot_2_files("test"+str(i))
