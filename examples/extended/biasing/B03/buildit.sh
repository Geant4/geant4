#!/bin/tcsh
source env.csh
cd include/
echo Calling swig
set SWIG=/opt/Anaphe/5.0.1/specific/redhat72/gcc-2.95.2/PublicDomainPackages/2.0.0/bin/swig-1.3.15
$SWIG -python -c++ -shadow -no_default B03App.i
echo moving files
mv B03App_wrap.cxx ../src/B03App_wrap.cc
mv B03App.py ../
cd ..
gmake

