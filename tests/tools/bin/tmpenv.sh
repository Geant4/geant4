if [ `uname -n | grep rsplus` ]; then
ulimit -d 400000
ulimit -s 100000
ulimit -m 100000
if [ X$REFTREE = Xprod ]; then
echo "Nothing set for TMPDIR!"
else
export TMPDIR=/afs/cern.ch/sw/geant4/stt/tmpdir
fi
fi

if [ `uname -n | grep sunasd1` ]; then
  echo "Nothing set for sunasd"
fi

if [ `uname -n | grep hpplus` ]; then
  echo "Nothing set for limit on hpplus"
# there is a pretty useless looking file from april in
# the directory below - I bet we are screwed if anything goes there
# with concurrent tests.
  if [ X$REFTREE = Xprod ]; then
    export TMPDIR=/afs/cern.ch/sw/geant4/stt/tmpdir
  else
    export TMPDIR=/afs/cern.ch/sw/geant4/stt/tmpdir
  fi
fi

if [ `uname -n | grep axcnsi` ]; then
  ulimit -d 200000
  ulimit -s 20000
  echo "axcnsi ulimit -d 200000 ulimit -s 20000"
  echo "Nothing set for TMPDIR on axcnsi!"
fi

if [ `uname -n | grep sgmedia` ]; then
  ulimit -s 200000
  ulimit -m 400000
  if [ X$REFTREE = Xprod ]; then
    export TMPDIR=/afs/cern.ch/sw/geant4/stt/tmpdir
  else
    export TMPDIR=/afs/cern.ch/sw/geant4/stt/tmpdir
  fi
fi
  
if [ `uname -n | grep dxplus` ]; then
  ulimit -d 400000
# (current understanding) TMPDIR is compiler loader  workspace (few 100 MBytes)
#                         G4TMP  are dependency files, object files  (up to Gbytes)
  if [ X$REFTREE = Xprod ]; then
    if [ ! -d /tmp/geant4tmp ]; then
      mkdir /tmp/geant4tmp
    fi
    export TMPDIR=/tmp/geant4tmp
#  else
#  export TMPDIR=/afs/cern.ch/sw/geant4/stt/tmpdir
  fi
fi
