#!/usr/local/bin/bash
##########################
# g4sbr.sh (SetupBuildRun)
##########################

#############set -x

REFTREE=$1
DEBOPT=`echo $2|cut -c 1`
NONISO=`echo $2 | grep NONISO`
REFTAG=$3
ACTION=$4
ACTARG1=$5
ACTARG2=$6
ACTARG3=$7
NONINCREMENTAL=$8

if [ $ACTION != all -a $ACTION != build -a $ACTION != run ]; then
  export G4LARGE_N=$ACTION
  ACTION=run
fi

if [ X$REFTREE = X -o X$DEBOPT = X -o X$REFTAG = X ]
then
  echo
  echo "Usage: g4sbr.sh dev[1|2]|prod debug[_ISO]|opt[_ISO] TAG_NAME Act arg1 arg2 arg3 NonIncrementalFlag"
  echo
  exit
fi
#

if [ $REFTREE = prod -o `echo $REFTREE | cut -c 1-3` = dev ]; then
:
else
  echo
  echo "Usage: First argument should be dev1/2 or prod."
  exit
fi

if [ X$DEBOPT = Xd -o X$DEBOPT = XD ]
then
  export G4DEBUG=1
elif [ X$DEBOPT = Xo -o X$DEBOPT = XO ]
then
  unset G4DEBUG
else
  echo
  echo "Usage: 2nd argument is debug or opt."
  exit
fi

if [ $NONISO ]
then
  export G4STTNONISO=1
else
  unset G4STTNONISO
fi

# Prevent compilation of histograms in examples code...
export G4NOHIST=1

####################################################################
# Setup environment in $REFTREE
####################################################################
echo "STT:SETUPEnvironment Complete"
# cd /afs/cern.ch/sw/geant4/stt/$REFTREE/testtools/geant4/tests/tools/bin
# so we can use pwd to determine the value of $REFTREE.
# . ${G4STTDIR}/bin/setup.sh (but this defines G4STTDIR).
cd /afs/cern.ch/sw/geant4/stt/$REFTREE/testtools/geant4/tests/tools/bin
.  /afs/cern.ch/sw/geant4/stt/$REFTREE/testtools/geant4/tests/tools/bin/setup.sh

env | grep G4
echo  "CLHEP_BASE_DIR $CLHEP_BASE_DIR"
ulimit -a
echo "STT:SETUPEnvironment Complete"

##########################
# Check if INPROGRESS
##########################
if [ -e $G4WORKDIR/inprogress.stat ]; then
  echo "STT:ABORT inprogress.stat exists."
  cat $G4WORKDIR/inprogress.stat
  ls -l $G4WORKDIR/inprogress.stat
  exit 
fi

###########################################
# Locks and stats
###########################################
trap "mv $G4WORKDIR/inprogress.stat $G4WORKDIR/interrupt.stat;touch $G4WORKDIR/interrupt.stat" TERM
rm $G4WORKDIR/done.stat $G4WORKDIR/interrupt.stat
touch $G4WORKDIR/inprogress.stat
cat >>  $G4WORKDIR/inprogress.stat <<EOF
${REFTAG}
EOF

######################################################################
# Prepare if not incremental: create new stt and link, clear bin|lib|tmp
######################################################################
echo "STT:SETUPDirectories Started"
if [ X$NONINCREMENTAL = X ]
then
  cd ${G4WORKDIR}/stt/${G4SYSTEM}
  NEXT_NUMBER=$[`ls -c1 gmake.log.*|sort|tail -1|cut -d "." -f3`+1]
  mv gmake.log gmake.log.${NEXT_NUMBER}
  echo "STT:UPDATE stt.${REFTAG} and RETAIN stt symbolic link."
else
  cd ${G4WORKDIR}
  if [ -d stt.${REFTAG} ]
  then
    echo "STT:ABORT stt.${REFTAG} already exists."
    exit
  fi
  echo "STT:CREATE stt.${REFTAG} and RESET stt symbolic link."
  mkdir stt.${REFTAG}
  rm -rf stt
  ln -s stt.${REFTAG} stt
  # To enable G4TMP to be in local/nfs disk and lib/bin in afs
  # (this is for public machines where geant4 does not own local disk space)
  # G4TMP is defined in specific.sh which is sourced in setup.sh
  G4TMP=${G4TMP:=$G4WORKDIR/tmp}
  export G4TMP
  G4LIB=${G4LIB:=$G4WORKDIR/lib}
  export G4LIB
  G4BIN=${G4BIN:=$G4WORKDIR/bin}
  export G4BIN
  echo "STT:REMOVE files in $G4TMP"
  rm -r $G4TMP/*
  echo "STT:REMOVE files in $G4LIB"
  rm -r $G4LIB/*
  echo "STT:REMOVE files in $G4BIN"
  rm -r $G4BIN/*
fi
echo "STT:SETUPDirectories Finished"
########################################################

################################
# Build&run all
################################
cd ${G4WORKDIR}
. ${G4STTDIR}/bin/limit.sh

if [ X$ACTION = Xbuild -o X$ACTION = Xall  ]
then
  . ${G4STTDIR}/bin/tmpenv.sh
  echo "STT:BUILD Started"
  ${G4STTDIR}/bin/build.sh $ACTARG1 $ACTARG2
  echo "STT:BUILD Finished"
fi

if [ X$ACTION = Xrun -o X$ACTION = Xall  ]
then
  echo "STT:RUN Started"
  ${G4STTDIR}/bin/run.sh $ACTARG3
  echo "STT:RUN Finished"
fi
####################################################################

########################
# Waiting for complete
########################
wait
mv $G4WORKDIR/inprogress.stat $G4WORKDIR/done.stat
touch $G4WORKDIR/done.stat

#######################################
# Mail about completion in TEST accont
#######################################
#mail serguei.sadilov@cern.ch <<EOF
#G4 sbr complete on $G4SYSTEM!
#EOF



