#!/usr/local/bin/bash
##########################
# g4sbr.sh (SetupBuildRun)
##########################

#############set -x

TREE=`echo $1|cut -c 1`
DEBOPT=`echo $2|cut -c 1`
REFTAG=$3
ACTION=$4
ACTARG1=$5
ACTARG2=$6
ACTARG3=$7
NONINCREMENTAL=$8
##############################################
# We have agreement about tag name: tag_name+
##############################################
PREVTAG=`echo $REFTAG|cut -d + -f1`

if [ X$TREE = X -o X$DEBOPT = X -o X$REFTAG = X ]
then
  echo
  echo "Usage: g4sbr.sh [dev|prod] [debug|opt] TAG_NAME Act arg1 arg2 arg3 NonIncrementalFlag"
  echo
  exit
fi
#

if [ X$TREE = Xd -o X$TREE = XD ]
then
  REFTREE=ref+
elif [ X$TREE = Xp -o X$TREE = XP ]
then
  REFTREE=ref
else
  echo
  echo "Usage: First argument is dev (uses ref+/) or prod (uses ref/)."
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

####################################################################
# Setup environment in $REFTREE
####################################################################
cd /afs/cern.ch/sw/geant4/stt/${REFTREE}/src/geant4/tests/tools/bin
. /afs/cern.ch/sw/geant4/stt/${REFTREE}/src/geant4/tests/tools/bin/setup.sh

env | grep G4
ulimit -a

##########################
# Check if INPROGRESS
##########################
if [ -e $G4WORKDIR/inprogress.stat ]; then
echo "In progress already!"
# exit 
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
# Prepare in ref[+] if not incremental: moving stt/, clear bin|lib|tmp
######################################################################
if [ X$NONINCREMENTAL = X ]
then
cd ${G4WORKDIR}/stt/${G4SYSTEM}
NEXT_NUMBER=$[`ls -c1 gmake.log.*|sort|tail -1|cut -d "." -f3`+1]
mv gmake.log gmake.log.${NEXT_NUMBER}
else
cd ${G4WORKDIR}
echo REMOVE
mv stt stt.${PREVTAG}
rm -r bin lib tmp
fi
########################################################

################################
# Build&run all in ref[+]
################################
cd ${G4WORKDIR}
. ${G4INSTALL}/tests/tools/bin/limit.sh

if [ X$ACTION = Xbuild -o X$ACTION = Xall  ]
then
#################
# Maybe workaround about first building sublibs without TMPDIR,
# but as afr as I know decision - use granular libs.
# Shortlived decision - first build ALL without TMPDIR - then
# with it.
################
#${G4INSTALL}/tests/tools/bin/build_specific.sh &
#${G4INSTALL}/tests/tools/bin/build.sh
. ${G4INSTALL}/tests/tools/bin/tmpenv.sh
${G4INSTALL}/tests/tools/bin/build.sh $ACTARG1 $ACTARG2
#unset $TMPDIR
#${G4INSTALL}/tests/tools/bin/build.sh test all
fi

#
# Shortlived solution for DEC6-AFS problems with templates
#
if [ X$G4SYSTEM = XDEC-cxx  ]
then
chmod +x ${G4WORKDIR}/bin/${G4SYSTEM}/*
fi

if [ X$ACTION = Xrun -o X$ACTION = Xall  ]
then
${G4INSTALL}/tests/tools/bin/run.sh $ACTARG3
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



