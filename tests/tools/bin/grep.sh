#
#  This scrip is used in conjunction with find to search
# something in the source code :
#
#  Examples :
#   UNIX> find $G4INSTALL/source -name "*h"  -exec ./grep.sh name {} \;
#   UNIX> find $G4INSTALL/source -name "*c"  -exec ./grep.sh name {} \;
#
grep $1 $2 && echo $2

