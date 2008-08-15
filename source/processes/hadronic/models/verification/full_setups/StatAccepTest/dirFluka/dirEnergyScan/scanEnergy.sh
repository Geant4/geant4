#------------------------------------------------------------------------
# Last update: 06-Aug-2008
#
# This simple shell script runs a energy scan, in the beam energy range
# 1 - 30 GeV, for 1 beam particle, and 4 types of simplified calorimeter
# (cms, atlas, pbwo4, tile).
#
# You have to pay attention to the creation of the executables.
# You need a different executable for each case, because they
# have different common block sizes, and also because you want
# to activate Birks in some cases.
# To build an executable, for instance for the CMS HCAL case,
# you have to do the following:
#    
#    $  cd dirUserRoutines
#    $  rm *.o 
#    $  rm mysmallinclude.inc
#    $  ln -s mysmallinclude.inc-cms mysmallinclude.inc
#    $  vi mgdraw.f    ->  decide to switch on/off Birks
#    $  cd ..
#    $  ./build
#    $  mv flukahp     flukahp-cms-birks
#    $  mv flukahp.map flukahp.map-cms-birks
#
# Notice, in particular, that you need to rename the 
# two files :  flukahp  and  flukahp.map .
# The renaming keeps the original filenames appending to
# them the same string  "-cms-birks" .
# The string that is appended to both  flukahp  and 
# flukahp  is specified by the following variable
# (that can be set below):
#   -$EXE_CMS    :  in the case of CMS HCAL
#   -$EXE_ATLAS  :  in the case of ATLAS HEC
#   -$EXE_PBWO4  :  in the case of PbWO4
#   -$EXE_TILE   :  in the case of TileCal
#
# Another variable is used in order to give a meaningful
# name to the final log file.
# The name is:
#          output.log-xxx-$PARTICLE-yyGeV-$LOG_zzz
# where:  "xxx" specified the type of calorimeters
#               (cms, atlas, pbwo4, tile);
#         "yy"  specified the beam energy in GeV
#               (1, 2, ... , 29, 30)
#         "zzz" specified the type of calorimeters
#               (CMS, ATLAS, PBWO4, TILE), where
#               LOG_zzz  is the last part of the
#               file name, to be specified below.
#
# See "***LOOKHERE***" below for the available choices.
#
#------------------------------------------------------------------------
#
#***LOOKHERE***
EXE_CMS=cms-birks ;    LOG_CMS=fluka-birks 
EXE_ATLAS=atlas ;      LOG_ATLAS=fluka
EXE_PBWO4=pbwo4 ;      LOG_PBWO4=fluka
EXE_TILE=tile-birks ;  LOG_TILE=fluka-birks
#
###PARTICLE=e
PARTICLE=pi
###PARTICLE=pi+
###PARTICLE=p
#**************
#
cd dirEnergyScan_$PARTICLE
#
echo " "
echo " === BEGIN === " FLUKA Energy Scan
echo " "
#
date
echo " === 1 GeV === "
cd dir1GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-1GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-1GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-1GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-1GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 2 GeV === "
cd dir2GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-2GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-2GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-2GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-2GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 3 GeV === "
cd dir3GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-3GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-3GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-3GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-3GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 4 GeV === "
cd dir4GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-4GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-4GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-4GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-4GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 5 GeV === "
cd dir5GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-5GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-5GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-5GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-5GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 6 GeV === "
cd dir6GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-6GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-6GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-6GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-6GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 7 GeV === "
cd dir7GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-7GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-7GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-7GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-7GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 8 GeV === "
cd dir8GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-8GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-8GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-8GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-8GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 9 GeV === "
cd dir9GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-9GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-9GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-9GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-9GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 10 GeV === "
cd dir10GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-10GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-10GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-10GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-10GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 11 GeV === "
cd dir11GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-11GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-11GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-11GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-11GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 12 GeV === "
cd dir12GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-12GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-12GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-12GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-12GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 13 GeV === "
cd dir13GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-13GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-13GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-13GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-13GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 14 GeV === "
cd dir14GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-14GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-14GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-14GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-14GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 15 GeV === "
cd dir15GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-15GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-15GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-15GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-15GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 16 GeV === "
cd dir16GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-16GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-16GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-16GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-16GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 17 GeV === "
cd dir17GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-17GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-17GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-17GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-17GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 18 GeV === "
cd dir18GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-18GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-18GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-18GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-18GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 19 GeV === "
cd dir19GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-19GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-19GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-19GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-19GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 20 GeV === "
cd dir20GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-20GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-20GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-20GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-20GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 21 GeV === "
cd dir21GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-21GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-21GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-21GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-21GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 22 GeV === "
cd dir22GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-22GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-22GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-22GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-22GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 23 GeV === "
cd dir23GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-23GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-23GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-23GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-23GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 24 GeV === "
cd dir24GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-24GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-24GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-24GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-24GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 25 GeV === "
cd dir25GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-25GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-25GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-25GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-25GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 26 GeV === "
cd dir26GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-26GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-26GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-26GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-26GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 27 GeV === "
cd dir27GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-27GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-27GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-27GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-27GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 28 GeV === "
cd dir28GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-28GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-28GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-28GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-28GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 29 GeV === "
cd dir29GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-29GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-29GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-29GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-29GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 30 GeV === "
cd dir30GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-30GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-30GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-30GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-30GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
echo " "
echo " === END === "
echo " "
