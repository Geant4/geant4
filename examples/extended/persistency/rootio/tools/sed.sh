#!/bin/sh
find . -name '*.h*' -exec ./g4exsed.sh "{}" \;
find . -name '*.h*' -exec ./g4psed.sh "{}" \;
find . -name '*.cc' -exec ./g4exsed.sh "{}" \;
find . -name '*.cc' -exec ./g4psed.sh "{}" \;
