#include "G4NeutronHPPhotonXSection.hh"

  G4NeutronHPPhotonXSection::G4NeutronHPPhotonXSection()
  {
    theExclusive = NULL;
    theExShell = NULL;
    theExEnergy = NULL;
    theExFlag = NULL;
    theExDisFlag = NULL;
  }
  G4NeutronHPPhotonXSection::~G4NeutronHPPhotonXSection()
  {
    if(theExclusive!=NULL) delete [] theExclusive;
    if(theExShell != NULL) delete [] theExShell;
    if(theExEnergy != NULL) delete [] theExEnergy;
    if(theExFlag != NULL) delete [] theExFlag;
    if(theExDisFlag != NULL) delete [] theExDisFlag;
  }
