#!/bin/sh
find . -name '*.hh' -exec $HOME/bin/g4pmv.sh "{}" \;
find . -name '*.cc' -exec $HOME/bin/g4pmv.sh "{}" \;
