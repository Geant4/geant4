#include "G4NeutronHPFinalState.hh"

  G4NeutronHPFinalState::G4NeutronHPFinalState()
  { 
    hasFSData = true; 
    hasXsec = true;
    hasAnyData = true;
    theBaseZ = 0;
    theBaseA = 0;
  };
  
  G4NeutronHPFinalState::~G4NeutronHPFinalState(){}
