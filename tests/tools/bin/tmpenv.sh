if [ `uname -n | grep rsplus` ]; then
ulimit -d 400000
ulimit -s 100000
ulimit -m 100000
if [ X$REFTREE = Xref ]; then
echo "Nothing set for TMPDIR!"
else
export TMPDIR=/afs/cern.ch/sw/geant4/stt/ref/Linux-g++/optim/other.tmp
fi
fi

if [ `uname -n | grep sunasd1` ]; then
echo "Nothing set!"
fi

if [ `uname -n | grep hpplus` ]; then
echo "Nothing set for limit!"
if [ X$REFTREE = Xref ]; then
export TMPDIR=/afs/cern.ch/sw/geant4/stt/ref/Linux-g++/optim/other.tmp
else
export TMPDIR=/afs/cern.ch/sw/geant4/stt/ref/Linux-g++/optim/other.tmp
fi
fi

if [ `uname -n | grep axcnsi` ]; then
ulimit -d 200000
ulimit -s 20000
echo "Nothing set for TMPDIR!"
fi

if [ `uname -n | grep sgmedia` ]; then
ulimit -s 200000
ulimit -m 400000
if [ X$REFTREE = Xref ]; then
export TMPDIR=/afs/cern.ch/sw/geant4/stt/ref/Linux-g++/optim/other.tmp
else
export TMPDIR=/afs/cern.ch/sw/geant4/stt/ref/Linux-g++/optim/other.tmp
fi
fi

if [ `uname -n | grep dxplus` ]; then
echo "I don't know about limits here :)...ooops...already know :("
ulimit -d 400000
fi
