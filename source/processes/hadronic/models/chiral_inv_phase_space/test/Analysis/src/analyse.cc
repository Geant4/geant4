#include "Analysis/src/ParticleInfo.h"
#include "globals.hh"
#include <fstream>

main()
{
  // ANAParticleInfo(G4double xSec, G4String aFileName);
  // Analyse the stuff
    // Inelastic cross-sections
    // p-Li6 = 128mb
    // p-Be  = 204mb
    // p-C   = 243mb
    // p-Al  = 455mb
    // p-Cu  = 823mb
    // p-Ta  = 1690mb
  
  ifstream steering("../inputs/analysis_steering");
  G4String fileName = "";
  while(steering >> fileName)
  {
    G4String path("../inputs/");
    G4String inputFile;
    inputFile = path+fileName;
    G4int flag;
    steering >> flag;
    G4cout << flag << G4endl;
    if(flag)
    {
      ifstream input(inputFile);
      G4cout << " IT! "<<inputFile<<G4endl;
      G4double crossSection;
      G4double A;
      G4double nevents;
      input >> crossSection;
      crossSection *= millibarn;
      input >> A;
      input >> nevents;
      G4double weight = crossSection/A;
    G4cout << "IT!!!!!"<<crossSection<<" "<<A<<" "<<nevents<<" "<<weight<< G4endl;
      G4String file("../logs/liste.");
      G4String it;
      it = file+fileName;
      ANAParticleInfo theInformation(weight, it);
      theInformation.Analyse();
      theInformation.Plot(fileName, nevents);
    }
    fileName = "";
  }
}
