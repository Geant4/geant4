#/bin/csh
#================================================
#     Macro file for Hadr05 run over all Physics Lists
#     26.06.2009 V.Ivanchneko
#================================================

mkdir -p $1
cd $1
rm -f *.*
setenv PHYSLIST $1
$G4BIN/$G4SYSTEM/Hadr05 $G4INSTALL/examples/extended/hadronic/Hadr05/pn.mac >& pn.out 
$G4BIN/$G4SYSTEM/Hadr05 $G4INSTALL/examples/extended/hadronic/Hadr05/pi.mac >& pi.out 
root -b -q  $G4INSTALL/examples/extended/hadronic/Hadr05/scripts/Plot.C
cd ../
#
