#include "G4NeutronHPAngularP.hh"

  G4NeutronHPAngularP::G4NeutronHPAngularP()
  {
    theCosTh = NULL;
    theProb = NULL;
  }
  G4NeutronHPAngularP::~G4NeutronHPAngularP()
  {
    if(theCosTh!=NULL) delete [] theCosTh;
    if(theProb!=NULL) delete [] theProb;
  }
  
