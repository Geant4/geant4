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
// $Id: G4EmSaturation.cc,v 1.6 2008-03-14 14:05:57 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmSaturation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.02.2008
//
// Modifications:
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmSaturation.hh"

#include "G4Proton.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmSaturation::G4EmSaturation()
:proton(0), tableManager(0), initialised(false)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmSaturation::~G4EmSaturation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmSaturation::Initialise()
{
 proton = G4Proton::Proton();
 tableManager = G4LossTableManager::Instance();
 initialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::BirksAttenuation(const G4Step* step) 
{
 if (!initialised) Initialise();
 //
 G4Material* material = step->GetTrack()->GetMaterial();
 G4double BirksConst  = material->GetIonisation()->GetBirksConstant();
 G4double destep      = step->GetTotalEnergyDeposit();
 
 if (BirksConst*destep <= 0.0) return destep;
 
 G4double nloss = step->GetNonIonizingEnergyDeposit();
 G4double eloss = std::max(destep - nloss, 0.);
 G4ParticleDefinition* particle = step->GetTrack()->GetDefinition(); 
 G4double z = std::abs(particle->GetPDGCharge()/eplus);

 // ionization 
 if (eloss*z > 0.) {
   G4double stepl = step->GetStepLength(); 
   eloss /= (1. + BirksConst*eloss/stepl);     
 }

 // non-ionizing loss 
 if (nloss*z > 0.) {
   G4double Anumber = ComputeAeff(material)/(g/mole);
   G4double escaled = nloss/Anumber;
   const
   G4MaterialCutsCouple* couple = step->GetTrack()->GetMaterialCutsCouple();
   G4double range = tableManager->GetRange(proton,escaled,couple)/(z*z); 
   nloss /= (1. + BirksConst*nloss/range);     
 }
 
 return (eloss + nloss);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmSaturation::ComputeAeff(G4Material* material)
{
  size_t nbElm                   = material->GetNumberOfElements();
  const G4ElementVector* element = material->GetElementVector();
  const G4double* massFract      = material->GetFractionVector();
  
  G4double Aeff = 0.0;
  for (size_t i=0; i<nbElm; i++) {
     Aeff += (massFract[i] * ((*element)[i]->GetA()));
  }
  
  return Aeff;   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
