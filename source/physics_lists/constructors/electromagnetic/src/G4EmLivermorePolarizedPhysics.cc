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

#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models

// gamma
#include "G4LivermorePolarizedPhotoElectricModel.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"

#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"

#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"

// interfaces
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"

// particles
#include "G4Gamma.hh"

//
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmLivermorePolarizedPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLivermorePolarizedPhysics::G4EmLivermorePolarizedPhysics(G4int ver, 
							     const G4String&)
  : G4EmLivermorePhysics(ver, "G4EmLivermorePolarized"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetEnablePolarisation(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLivermorePolarizedPhysics::~G4EmLivermorePolarizedPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLivermorePolarizedPhysics::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmLivermorePhysics::ConstructProcess();

  G4LossTableManager* man = G4LossTableManager::Instance();
  G4EmConfigurator* em_config = man->EmConfigurator();
  G4double livEnergyLimit = 1*GeV;

  // Add Livermore EM Processes
  G4VEmModel* mod = new G4LivermorePhotoElectricModel();
  mod->SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized());
  em_config->SetExtraEmModel("gamma", "phot", mod);

  G4VEmModel* comptLiv = new G4LivermorePolarizedComptonModel();
  comptLiv->SetHighEnergyLimit(livEnergyLimit);
  em_config->SetExtraEmModel("gamma", "compt", comptLiv);

  G4VEmModel* convLiv = new G4LivermorePolarizedGammaConversionModel();
  convLiv->SetHighEnergyLimit(livEnergyLimit);
  em_config->SetExtraEmModel("gamma", "conv", convLiv);

  G4VEmModel* theRay = new G4LivermorePolarizedRayleighModel();
  em_config->SetExtraEmModel("gamma", "Rayl", theRay);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
