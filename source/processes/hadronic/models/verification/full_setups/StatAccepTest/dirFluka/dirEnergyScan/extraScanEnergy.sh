#------------------------------------------------------------------------
# Last update: 21-Aug-2008
#
# This simple shell script runs some beam energy points
# above 30 GeV, for 1 beam particle, and 4 types of 
# simplified calorimeter (cms, atlas, pbwo4, tile).
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
#
date
echo " === 50 GeV === "
cd dir50GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-50GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-50GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-50GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-50GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
#
date
echo " === 100 GeV === "
cd dir100GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-100GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-100GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-100GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-100GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
#
###date
###echo " === 180 GeV === "
###cd dir180GeV
###rm *00*
####
###ln -sfn ../../flukahp-$EXE_CMS     flukahp
###ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
###mv cmsHCAL001.log output.log-cms-$PARTICLE-180GeV-$LOG_CMS 
####
###ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
###ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
###mv atlasHEC001.log output.log-atlas-$PARTICLE-180GeV-$LOG_ATLAS
####
###ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
###ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
###mv pbwo4001.log output.log-pbwo4-$PARTICLE-180GeV-$LOG_PBWO4
####
###ln -sfn ../../flukahp-$EXE_TILE     flukahp
###ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
###mv tileCal001.log output.log-tile-$PARTICLE-180GeV-$LOG_TILE
####
###rm flukahp flukahp.map
###cd ..
#
#
date
echo " === 200 GeV === "
cd dir200GeV
rm *00*
#
ln -sfn ../../flukahp-$EXE_CMS     flukahp
ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
mv cmsHCAL001.log output.log-cms-$PARTICLE-200GeV-$LOG_CMS 
#
ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
mv atlasHEC001.log output.log-atlas-$PARTICLE-200GeV-$LOG_ATLAS
#
ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
mv pbwo4001.log output.log-pbwo4-$PARTICLE-200GeV-$LOG_PBWO4
#
ln -sfn ../../flukahp-$EXE_TILE     flukahp
ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
mv tileCal001.log output.log-tile-$PARTICLE-200GeV-$LOG_TILE
#
rm flukahp flukahp.map
cd ..
#
#
###date
###echo " === 300 GeV === "
###cd dir300GeV
###rm *00*
####
###ln -sfn ../../flukahp-$EXE_CMS     flukahp
###ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
###mv cmsHCAL001.log output.log-cms-$PARTICLE-300GeV-$LOG_CMS 
####
###ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
###ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
###mv atlasHEC001.log output.log-atlas-$PARTICLE-300GeV-$LOG_ATLAS
####
###ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
###ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
###mv pbwo4001.log output.log-pbwo4-$PARTICLE-300GeV-$LOG_PBWO4
####
###ln -sfn ../../flukahp-$EXE_TILE     flukahp
###ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
###mv tileCal001.log output.log-tile-$PARTICLE-300GeV-$LOG_TILE
####
###rm flukahp flukahp.map
###cd ..
#
#
###date
###echo " === 350 GeV === "
###cd dir350GeV
###rm *00*
####
###ln -sfn ../../flukahp-$EXE_CMS     flukahp
###ln -sfn ../../flukahp.map-$EXE_CMS flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp cmsHCAL
###mv cmsHCAL001.log output.log-cms-$PARTICLE-350GeV-$LOG_CMS 
####
###ln -sfn ../../flukahp-$EXE_ATLAS     flukahp
###ln -sfn ../../flukahp.map-$EXE_ATLAS flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp atlasHEC
###mv atlasHEC001.log output.log-atlas-$PARTICLE-350GeV-$LOG_ATLAS
####
###ln -sfn ../../flukahp-$EXE_PBWO4     flukahp
###ln -sfn ../../flukahp.map-$EXE_PBWO4 flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp pbwo4
###mv pbwo4001.log output.log-pbwo4-$PARTICLE-350GeV-$LOG_PBWO4
####
###ln -sfn ../../flukahp-$EXE_TILE     flukahp
###ln -sfn ../../flukahp.map-$EXE_TILE flukahp.map
###$FLUPRO/flutil/rfluka -M1 -N0 -e flukahp tileCal
###mv tileCal001.log output.log-tile-$PARTICLE-350GeV-$LOG_TILE
####
###rm flukahp flukahp.map
###cd ..
#
#
date
#
echo " "
echo " === END === "
echo " "
