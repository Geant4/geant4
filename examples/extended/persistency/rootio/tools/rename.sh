#!/bin/sh
find . -name '*.h*' -exec ./g4exmv.sh "{}" \;
find . -name '*.cc' -exec ./g4exmv.sh "{}" \;
