#include "G4NeutronHPLegendreTable.hh"

  G4NeutronHPLegendreTable::G4NeutronHPLegendreTable()
  {
    nCoeff=0; 
    theCoeff = NULL;
  }
  G4NeutronHPLegendreTable::~G4NeutronHPLegendreTable()
  {if(theCoeff!=NULL) delete [] theCoeff;}
