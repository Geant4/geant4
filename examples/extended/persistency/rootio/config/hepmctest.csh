#! /bin/csh -f

mkdir -p $G4PEX_DIR/tmp/$G4SYSTEM/.hepmc
set test_src = $G4PEX_DIR/tmp/$G4SYSTEM/.hepmc/qhepmc.cc
set test_obj = $G4PEX_DIR/tmp/$G4SYSTEM/.hepmc/qhepmc.o

cat > $test_src  << +EOF
#include "HepMC/GenVertex.h"

int main()
{  
  HepMC::GenVertex av;
  av.barcode();
}
+EOF

g++ -c -I$HEPMC_DIR -I$CLHEP_BASE_DIR/include $test_src \
    -o $test_obj >& /dev/null
set qbar=$status

\rm -r $G4PEX_DIR/tmp/$G4SYSTEM/.hepmc

if ( $qbar == 0 ) then
  echo true
else
  echo false
endif

