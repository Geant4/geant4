#!/bin/sh

class="$1"
dict="$2"

if [ "$class" = "" ]; then
  echo "Usage: `basename $0`  'class names' [ dict_file ]"
  exit 1
fi

if [ "$dict" = "" ]; then
  dict="DefaultLinkDef.h"
fi

echo "Creating ROOTCINT LinkDef $2 for classes $1"

echo "#ifdef __CINT__"                  > $dict
echo ""                                >> $dict
echo "#pragma link off all globals;"   >> $dict
echo "#pragma link off all classes;"   >> $dict
echo "#pragma link off all functions;" >> $dict
echo ""                                >> $dict

for name in $class ; do 
 echo "#pragma link C++ class ${name};" >> $dict
done

echo ""                                >> $dict
echo "#endif"                          >> $dict

