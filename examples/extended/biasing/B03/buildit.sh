#!/bin/tcsh
# The script uses swig to create c++ code according to the 
# file include/B03App.i. It also creates the file B03App.py 
# containing the python class wrappers. 
# If you need to use this script you have to set the 
# variable SWIG according to your installation of swig.

# If you want to use this script please edit the file
# env.csh according to your installation of Anaphe and source
# it before running this script.

source env.csh
cd include/
echo Calling swig
set SWIG=${ANAPHE_BASE}/PublicDomainPackages/2.0.0/bin/swig-1.3.15
$SWIG -python -c++ -shadow -no_default B03App.i
echo moving files
mv B03App_wrap.cxx ../src/B03App_wrap.cc
mv B03App.py ../
cd ..
gmake

