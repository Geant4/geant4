####################################################
# Works for bash, ksh and zsh
# Usage :
#    [ba,k,z]sh> cd <g4install>/tests/tools/bin
#    [ba,k,z]sh> . setup.sh
####################################################
#

#########set -x

# Execute stt members specifc setup first :
found=no
if [ `pwd | grep /stt/dev1/` ]; then
  . /afs/cern.ch/sw/geant4/stt/dev1/src/geant4/tests/tools/bin/specific.sh
  found=yes
fi
if [ `pwd | grep /stt/dev2/` ]; then
  . /afs/cern.ch/sw/geant4/stt/dev2/src/geant4/tests/tools/bin/specific.sh
  found=yes
fi
if [ `pwd | grep /stt/prod/` ]; then
  . /afs/cern.ch/sw/geant4/stt/prod/src/geant4/tests/tools/bin/specific.sh
  found=yes
fi
if [ $found = no ]; then
  echo Cannot find sensible tree - pwd should contain /stt/prod/, etc.
  exit
fi

if [ -z "$G4SYSTEM" ]; then
  echo "You have first to set environment variable G4SYSTEM !"
  exit
else
  echo "You are working on a $G4SYSTEM system"
#...

# Some checks :
if [ -z "$G4INSTALL" ]; then
 echo "You have first to set environment variable G4INSTALL !"
 exit
else
 echo "You use for the source $G4INSTALL"
#...

if [ -z "$G4WORKDIR" ]; then
 echo "You have first to set environment variable G4WORKDIR !"
 exit
else
 echo "You use as a work directory $G4WORKDIR "
fi
#
if [ -z "$G4LIB" ]; then
 echo "You have first to set environment variable G4LIB !"
 exit
else
 echo "You are using as a library directory $G4LIB "
fi
#
# Other G4 environment variables.
#
export NeutronHPCrossSections=$G4INSTALL/../G4NDL0.2
export G4LEVELGAMMADATA=$G4INSTALL/data/PhotonEvaporation
export G4LEDATA=$G4INSTALL/../G4EMLOW0.3
#
# Some alias -kes :
alias g4root="cd $G4INSTALL"
alias g4source="cd $G4INSTALL/source"
alias g4config="cd $G4INSTALL/config"
alias g4vis="cd $G4INSTALL/source/visualization"
alias g4tmp="cd $G4WORKDIR/tmp/$G4SYSTEM"
alias g4bin="cd $G4WORKDIR/bin/$G4SYSTEM"
alias g4lib="cd $G4LIB/$G4SYSTEM"
alias g4stt="cd $G4WORKDIR/stt/$G4SYSTEM"
#
alias g4tests="cd $G4INSTALL/tests"
alias g4tools="cd $G4INSTALL/tests/tools/bin"
alias g4nt="cd $G4INSTALL/tests/tools/NT"
#
# Below aliases assume that $G4WORKDIR/stt/$G4SYSTEM exists !
alias g4make="gmake global> & $G4WORKDIR/stt/$G4SYSTEM/gmake.log &"
alias g4build="$G4INSTALL/tests/tools/bin/build.sh"
alias g4run="$G4INSTALL/tests/tools/bin/run.sh"
alias g4analyse="$G4INSTALL/tests/tools/bin/analyse.sh"
alias g4tail="tail -f $G4WORKDIR/stt/$G4SYSTEM/gmake.log"
alias g4filter="$G4INSTALL/tests/tools/bin/filter.sh $G4SYSTEM | more"
#
#alias search="$G4INSTALL/tests/search.sh"
# Make sure the test tool lib exists
if [ ! -d $G4INSTALL/tests/tools/lib/$G4SYSTEM ];then
 mkdir $G4INSTALL/tests/tools/lib/$G4SYSTEM
 cd  $G4INSTALL/tests/tools/lib; gmake
 cd  $G4INSTALL/tests/tools/bin
else
 if [ -n `ls $G4INSTALL/tests/tools/lib/$G4SYSTEM/lib* ` ];then 
  cd  $G4INSTALL/tests/tools/lib; gmake
  cd  $G4INSTALL/tests/tools/bin
 fi
fi

fi # if [ -z "$G4INSTALL" ]

fi # if [ -z "$G4SYSTEM" ]