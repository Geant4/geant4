####################################################
#
# Usage :
#     csh> cd <g4install>/tests/tools/bin
#     csh> source setup.csh
####################################################
#
# Execute stt members specifc setup first :
source specific.csh
#
# Some checks :
if ( $?G4INSTALL ) then
else
  echo "You have first to set environment variable G4INSTALL !"
  exit
endif
#
if ( $?G4SYSTEM ) then
else
  echo "You have first to set environment variable G4SYSTEM !"
  exit
endif
#
if ( $?G4WORKDIR ) then
else
  echo "You have first to set environment variable G4WORKDIR !"
  exit
endif
#
if ( $?G4LIB ) then
else
  echo "You have first to set environment variable G4LIB !"
  exit
endif
#
# Other G4 environment variables.
#
setenv NeutronHPCrossSections $G4INSTALL/../G4NDL0.2
setenv G4LEVELGAMMADATA $G4INSTALL/data/PhotonEvaporation
setenv G4RADIOACTIVEDATA $G4INSTALL/data/RadiativeDecay
setenv G4LEDATA $G4INSTALL/../G4EMLOW0.3
#
# Some aliases :
alias g4root   "cd $G4INSTALL"
alias g4source "cd $G4INSTALL/source"
alias g4config "cd $G4INSTALL/config"
alias g4vis    "cd $G4INSTALL/source/visualization"
alias g4tmp    "cd $G4WORKDIR/tmp/$G4SYSTEM"
alias g4bin    "cd $G4WORKDIR/bin/$G4SYSTEM"
alias g4lib    "cd $G4LIB/$G4SYSTEM"
alias g4stt    "cd $G4WORKDIR/stt/$G4SYSTEM"
#
alias g4tests   "cd $G4INSTALL/tests"
alias g4test201 "cd $G4INSTALL/tests/test201"
alias g4tools   "cd $G4INSTALL/tests/tools/bin"
alias g4nt      "cd $G4INSTALL/tests/tools/NT"
alias g4n02     "cd $G4INSTALL/examples/novice/N02"
alias g4n03     "cd $G4INSTALL/examples/novice/N03"
#
# Below aliases assume that $G4WORKDIR/stt/$G4SYSTEM exists !
alias g4make    "gmake global> & $G4WORKDIR/stt/$G4SYSTEM/gmake.log &"
alias g4build   "$G4INSTALL/tests/tools/bin/build.sh"
alias g4run     "$G4INSTALL/tests/tools/bin/run.sh"
alias g4analyse "$G4INSTALL/tests/tools/bin/analyse.sh"
alias g4tail    "tail -f $G4WORKDIR/stt/$G4SYSTEM/gmake.log"
alias g4filter  "$G4INSTALL/tests/tools/bin/filter.sh $G4SYSTEM | more"
#
alias test201 "$G4WORKDIR/bin/$G4SYSTEM/test201"
#
#alias search "$G4INSTALL/tests/search.sh"
