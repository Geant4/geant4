#include "G4NeutronHPPolynomExpansion.hh"

  G4NeutronHPPolynomExpansion::G4NeutronHPPolynomExpansion()
  {
    theCoeff = NULL;
    nPoly=0;
  }
  
  G4NeutronHPPolynomExpansion::~G4NeutronHPPolynomExpansion()
  {
    if(theCoeff!=NULL) delete [] theCoeff;
  }
