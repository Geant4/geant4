#!/bin/tcsh
source env.csh
cd include/
echo Calling swig
swig -python -c++ -shadow -no_default B03App.i
echo moving files
mv B03App_wrap.cxx ../src/B03App_wrap.cc
mv B03App.py ../
cd ..
gmake

