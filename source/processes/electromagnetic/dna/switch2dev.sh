#!/bin/bash

CMF=`find . -name 'CMakeLists.txt'`


for file in $CMF ; do
svn merge $file -r 80081:80074
done
