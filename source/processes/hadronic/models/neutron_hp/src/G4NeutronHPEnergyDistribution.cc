#include "G4NeutronHPEnergyDistribution.hh"

  G4NeutronHPEnergyDistribution::G4NeutronHPEnergyDistribution()
  {
    theEnergyDistribution = NULL;
    theNumberOfPartials = 0;
    theRepresentationType = 0;
  }
  G4NeutronHPEnergyDistribution::~G4NeutronHPEnergyDistribution()
  {
    if(theEnergyDistribution != NULL)
    {
      for(G4int i=0; i<theNumberOfPartials; i++) 
      {
        delete theEnergyDistribution[i];
      }
      delete [] theEnergyDistribution;
    }
  }
