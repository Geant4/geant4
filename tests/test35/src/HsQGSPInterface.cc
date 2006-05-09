//----------------------------------------------------------------------------
//
//  Package   : Simulation 
//
// Description: Algorithm of G4 HARP for Hadron Production in the target
//
// Author:      V.Ivanchenko 05.03.04
//
// Modifications: 
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HsQGSPInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

HsQGSPInterface::HsQGSPInterface() 
{
  theModel = new G4TheoFSGenerator;
  theCascade = new G4GeneratorPrecompoundInterface;
  thePreEquilib = new G4PreCompoundModel(&theHandler);
  theCascade->SetDeExcitation(thePreEquilib);  
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(&theStringModel);
  theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
  theModel->SetMinEnergy(6.*GeV);
  theModel->SetMaxEnergy(100*TeV);
}

HsQGSPInterface::~HsQGSPInterface() 
{
  delete theStringDecay;
}

G4HadFinalState* HsQGSPInterface::ApplyYourself(const G4HadProjectile& aTrack, 
                                                        G4Nucleus& theNucleus)
{
  return theModel->ApplyYourself(aTrack, theNucleus);
}


