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
//
// GEANT4 Class file
//
//
// File name:     G4SimplePsAtRestModel
//
// Author:        I.Semeniouk & D.Bernard
//
// Creation date: 04 Juin 2024
//
// -------------------------------------------------------------------
//

#include "G4SimplePsAtRestModel.hh"
#include "G4SimplePositronAtRestModel.hh"
#include "G4OrePowellAtRestModel.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SimplePsAtRestModel::G4SimplePsAtRestModel()
  : G4VPositronAtRestModel("SimplePs")
{
  f3gFranction = G4EmParameters::Instance()->OrtoPsFraction();
  model2g = new G4SimplePositronAtRestModel();
  model3g = new G4OrePowellAtRestModel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SimplePsAtRestModel::SampleSecondaries(
             std::vector<G4DynamicParticle*>& secParticles,
             G4double& localEnergyDeposit, const G4Material* mat) const
{
  // G4cout << "SampleSecondaries model " << GetName() << G4endl;
  // G4cout << "3 gamma fraction " << f3gFranction  << G4endl;
  if ( G4UniformRand() >  f3gFranction ) {
    model2g->SampleSecondaries(secParticles,localEnergyDeposit,mat);
  } else {
    model3g->SampleSecondaries(secParticles,localEnergyDeposit,mat);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SimplePsAtRestModel::PrintGeneratorInformation() const
{
  G4cout << G4endl;
  model2g->PrintGeneratorInformation();
  model3g->PrintGeneratorInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SimplePsAtRestModel::~G4SimplePsAtRestModel()
{
  delete model2g;
  delete model3g;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
