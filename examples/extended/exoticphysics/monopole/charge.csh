#/bin/csh
# 
# 23 March 2010 V.Ivanchenko test on heavy objects with great charge
#

gmake

mv res+.out res+.out.old
mv res-.out res-.out.old

$G4BIN/Linux-g++/monopole monopole.in '0 100 150 GeV'  >& res+.out
$G4BIN/Linux-g++/monopole monopole.in '0 -100 150 GeV' >& res-.out

tkdiff res+.out res-.out

#
