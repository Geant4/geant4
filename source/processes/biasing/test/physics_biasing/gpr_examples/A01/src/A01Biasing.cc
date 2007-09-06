#include "A01Biasing.hh"

//#include "G4GPRBuilder.hh"
#include "G4GPRConverter.hh"
#include "G4Electron.hh"

void A01Biasing::ConstructBiasing() 
{
  G4cout<<"jane construct a01 biasing"<<G4endl;
  //G4GPRConverter::LoadGPR();
  //  G4GPRConverter::LoadGPR(G4Electron::Definition());
  
  // Override default A01 physics list for positrons. Just for fun, make the positrons use
  // the low energy physics list by default
  G4String newPhysicsList("New Default");
  CreateDefaultPhysicsList<G4Positron>(newPhysicsList);
  
  AddProcess<G4Positron>(new G4Transportation, x, newPhysicsList);
  AddProcess<G4Positron>(new G4MultipleScattering, -1, 1, 1, newPhysicsList);
  AddProcess<G4Positron>(new G4eIonisation,        -1, 2, 2, newPhysicsList);
  AddProcess<G4Positron>(new G4eBremsstrahlung,    -1, 3, 3, newPhysicsList);
  AddProcess<G4Positron>(new G4eplusAnnihilation,   0,-1, 4, newPhysicsList);
}
