#------------------------------------------------------------------------
# Last update: 04-Sep-2008
#
# This simple shell script runs a energy scan, in the beam energy range
# 1 - 30 GeV, and beyond, for 1 beam particle, and 1 type of 
# simplified calorimeter, between the following 6 choices:
#  cms , atlas , pbwo4 , tile , atile , precisecms .
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
# flukahp  is specified by the $EXE variable (see below).
#
# Another variable is used in order to give a meaningful
# name to the final log file.
# The name is:
#          output.log-${TYPE}-${PARTICLE}-yyGeV-$LOG
# where:  $TYPE specified the type of calorimeters
#               (cms, atlas, pbwo4, tile, atile, precisecms);
#         "yy"  specified the beam energy in GeV
#               (1, 2, ... , 29, 30, ...)
#         $LOG  is the last part of the file name, 
#               to be specified below.
#
# See "***LOOKHERE***" below for the available choices.
#
# This script is similar to scanEnergy1.sh, with two important
# differences:
#  a) it allows to run only one type of calorimeters, between
#     6 possibilities 
#     (whereas scanEnergy.sh runs 4 types of calorimeters: 
#      cms, atlas, pbwo4, tile);
#  b) it runs the all spectrum of beam energies
#     (whereas scanEnergy.sh runs between 1 and 30 GeV, and
#      then you need extraScanEnergy.sh to runs beyond that,
#      always for 4 calorimeters: cms, atlas, pbwo4, tile).
#
#------------------------------------------------------------------------
#
#***LOOKHERE***
###TYPE=cms
###EXE=$TYPE-birks ;  LOG=fluka-birks
###NAME_INP=cmsHCAL
#
###TYPE=atlas
###EXE=$TYPE ;  LOG=fluka
###NAME_INP=atlasHEC
#
###TYPE=pbwo4
###EXE=$TYPE ;  LOG=fluka
###NAME_INP=pbwo4
#
###TYPE=tile
###EXE=$TYPE-birks ;  LOG=fluka-birks
###NAME_INP=tileCal
#
TYPE=atile
EXE=$TYPE-birks ;  LOG=fluka-birks
NAME_INP=atlasTILE
#
###TYPE=precisecms
###EXE=$TYPE-birks ;  LOG=fluka-birks
###NAME_INP=precise_cmsHCAL
#
###PARTICLE=e
PARTICLE=pi
###PARTICLE=pi+
###PARTICLE=p
#**************
#
cd dirEnergyScan_${PARTICLE}
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
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-1GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 2 GeV === "
cd dir2GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-2GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 3 GeV === "
cd dir3GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-3GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 4 GeV === "
cd dir4GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-4GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 5 GeV === "
cd dir5GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-5GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 6 GeV === "
cd dir6GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-6GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 7 GeV === "
cd dir7GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-7GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 8 GeV === "
cd dir8GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-8GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 9 GeV === "
cd dir9GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-9GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 10 GeV === "
cd dir10GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-10GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 11 GeV === "
cd dir11GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-11GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 12 GeV === "
cd dir12GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-12GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 13 GeV === "
cd dir13GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-13GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 14 GeV === "
cd dir14GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-14GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 15 GeV === "
cd dir15GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-15GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 16 GeV === "
cd dir16GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-16GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 17 GeV === "
cd dir17GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-17GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 18 GeV === "
cd dir18GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-18GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 19 GeV === "
cd dir19GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-19GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 20 GeV === "
cd dir20GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-20GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 21 GeV === "
cd dir21GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-21GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 22 GeV === "
cd dir22GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-22GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 23 GeV === "
cd dir23GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-23GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 24 GeV === "
cd dir24GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-24GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 25 GeV === "
cd dir25GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-25GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 26 GeV === "
cd dir26GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-26GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 27 GeV === "
cd dir27GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-27GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 28 GeV === "
cd dir28GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-28GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 29 GeV === "
cd dir29GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-29GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
date
echo " === 30 GeV === "
cd dir30GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-30GeV-${LOG}
#
rm flukahp flukahp.map
cd ..
#
#
date
echo " === 50 GeV === "
cd dir50GeV
rm *00*
#
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-50GeV-${LOG}
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
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-100GeV-${LOG}
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
###ln -sfn ../../flukahp-${EXE}     flukahp
###ln -sfn ../../flukahp.map-${EXE} flukahp.map
###${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
###mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-180GeV-${LOG}
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
ln -sfn ../../flukahp-${EXE}     flukahp
ln -sfn ../../flukahp.map-${EXE} flukahp.map
${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-200GeV-${LOG}
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
###ln -sfn ../../flukahp-${EXE}     flukahp
###ln -sfn ../../flukahp.map-${EXE} flukahp.map
###${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
###mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-300GeV-${LOG}
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
###ln -sfn ../../flukahp-${EXE}     flukahp
###ln -sfn ../../flukahp.map-${EXE} flukahp.map
###${FLUPRO}/flutil/rfluka -M1 -N0 -e flukahp ${NAME_INP}
###mv ${NAME_INP}001.log output.log-${TYPE}-${PARTICLE}-350GeV-${LOG}
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
