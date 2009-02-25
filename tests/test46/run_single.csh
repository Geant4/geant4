#/bin/csh

rm -f $1.out

echo 'Run test46 at ' $HOST >& $1.out

$G4MY/test46 $G4INSTALL/tests/test46/$1.in >>& $1.out
mv test46.root $1.root

#
