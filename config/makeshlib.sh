#!/bin/sh -f
############################################
##### To build GEANT4 shared libraries #####
############# G.Barrand ####################
############################################
if test $# -eq 2 ; then
g4system=$1
name=$2
liba=lib$name.a
#echo $g4system
#echo $name
############################################
##### build libs ###########################
############################################
if test -f $g4system/$liba ; then

if test `uname` = "OSF1" ; then
cd $g4system
#  Need to specify a path, so that at run time loader
# is able to find the so lib. This path is put in the 
# binaries at link time. In principle a relative
# path from the binary is ok, but we will have
# more than one binary at different places; so
# we have to specify an absolute path.
rpath=`pwd`
#echo $rpath
libso=lib$name.so
/bin/rm -f $libso
/bin/rm -f so_locations
ld -expect_unresolved '*' -shared -o $libso -all $liba -rpath $rpath
cd ..
echo "$libso built."
fi

if test `uname` = "HP-UX" ; then
cd $g4system
libso=lib$name.sl
/bin/rm -f $libso
# get .o files
objs=`ls *.o`
#ld -b -a shared -o $libso $objs
# CC does more than ld.
CC -b -o $libso $objs
cd ..
echo "$libso built."
fi

if test `uname` = "IRIX" ; then
cd $g4system
# For path same remarks as for OSF1.
rpath=`pwd`
libso=lib$name.so
/bin/rm -f $libso
/bin/rm -f so_locations
#ld -shared -o $libso -all $liba  -rpath $rpath
# CC does more than ld.
CC -Wl,-rpath,$rpath -shared -o $libso -all $liba
cd ..
echo "$libso built."
fi

if test `uname` = "Linux" ; then
cd $g4system
# For path same remarks as for OSF1.
libso=lib$name.so
/bin/rm -f $libso
g++ -Wl,-soname,$libso -shared *.o -o $libso
echo "$libso built."
fi

if test `uname` = "AIX" ; then
# Not yet ready.
exit
#set -x
cd $g4system
libso=lib$name'.so'
/bin/rm -f lib$name.o
/bin/rm -f $libso
/bin/rm -f lib$name.exp
# Get .o files
objs=`ls *.o`
# Determine how to invoke nm depending on AIX version
AIXVERSION=`uname -v`
case ${AIXVERSION} in
#{
        3*)
                NM=/usr/ucb/nm
                ;;
        4*)
                NM='/usr/bin/nm -B'
                ;;
        *)
                echo "Error in g4makeshlib.sh!"
                exit 1
                ;;
#}
esac
# Build export file (from Mesa mklib.aix).
echo "#! lib$name.so " > lib$name.exp
echo "noentry" >> lib$name.exp
$NM $objs | awk '/ [BD] /{print $3}' | sort | uniq >> lib$name.exp
# Below command from Mesa mklib.aix.
#cc -o lib$name.o $objs $3 -bE:lib$name.exp -bM:SRE -enoentry 
# Unnecessary option -enoentry that produces
#  noentry file that contains exported symbols
cc -o lib$name.o $objs $3 -bE:lib$name.exp -bM:SRE
if test -f lib$name.o ; then
ar r $libso lib$name.o
/bin/rm -f lib$name.o
else
echo "Can't create lib$name.o."
fi
#cp  $libso lib$name.a
cd ..
echo "$libso treated."
fi


fi
############################################
else
echo 'Give name of binaray directory and library'
fi
