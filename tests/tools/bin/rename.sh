#!/bin/sh

if [ $# -lt 3 ]
then
   echo Usage: rename.sh oldstring newstring files
   exit 1
fi

echo Replacing 1st occurrence of $1 with $2 \
     in the following files:

a=$1; shift;
b=$1; shift;
for i
do
   j=`echo $i | sed s/$a/$b/`
   if [ $i = $j ]
   then
      echo Nothing to replace in $i
   else
      echo -n "Renaming $i to $j - OK? (y/n): "
      read OK
      case $OK in
         y*|Y*)  mv $i $j
                 echo $i renamed to $j;;
         *)      ;;
      esac
   fi
done
exit 0
