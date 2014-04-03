#!/bin/bash

CMF=`find . -name 'CMakeLists.txt'`


for file in $CMF ; do
svn merge $file -r 80074:80081
done
