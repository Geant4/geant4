//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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
  theModel->SetMinEnergy(GeV);
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


