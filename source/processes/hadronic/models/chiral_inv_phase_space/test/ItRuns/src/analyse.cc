#include "Analysis/src/ParticleInfo.h"
#include "globals.hh"

main()
{
  // ANAParticleInfo(G4double xSec, G4String aFileName);
  // Analyse the stuff
  ANAParticleInfo theInformation(100.*millibarn, "../logs/liste");
  theInformation.Analyse();
}
