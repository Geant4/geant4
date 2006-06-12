//----------------------------------------------------------------------------
//
//  Package   : Simulation 
//
// Description: Algorithm of G4 HARP for Hadron Production in the target
//
// Author:      V.Ivanchenko 12.06.06
//
// Modifications: 
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HsFTFCInterface.hh"
#include "G4LundStringFragmentation.hh"
#include "G4FTFModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

HsFTFCInterface::HsFTFCInterface() 
{
  theModel = new G4TheoFSGenerator;
  theStringModel = new G4FTFModel();
  theCascade = new G4StringChipsParticleLevelInterface;
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(theStringModel);
  theStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation());
  theStringModel->SetFragmentationModel(theStringDecay);
  theModel->SetMinEnergy(6.*GeV);
  theModel->SetMaxEnergy(100*TeV);
}

HsFTFCInterface::~HsFTFCInterface() 
{
  delete theModel;
}

G4HadFinalState* HsFTFCInterface::ApplyYourself(const G4HadProjectile& aTrack, 
                                                        G4Nucleus& theNucleus)
{
  return theModel->ApplyYourself(aTrack, theNucleus);
}


