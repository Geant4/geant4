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
// $Id: G4ePolarizedBremsstrahlung.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ePolarizedBremsstrahlung
//
// Author:        Karim Laihem
//
// Creation date: 26.06.2005
//
// Modifications:
//    19-08-06 addapted to accomodate geant481 structure
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ePolarizedBremsstrahlung.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4ePolarizedBremsstrahlungModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4ePolarizedBremsstrahlung::G4ePolarizedBremsstrahlung(const G4String& name):
  G4eBremsstrahlung(name)
{}

void G4ePolarizedBremsstrahlung::InitialiseEnergyLossProcess(
            const G4ParticleDefinition*,
            const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    SetSecondaryParticle(G4Gamma::Gamma());
    SetIonisation(false);

    G4VEmFluctuationModel* fm =  nullptr;
    //G4VEmFluctuationModel* fm = new G4UniversalFluctuation();

    G4VEmModel* em = new G4ePolarizedBremsstrahlungModel;
    em->SetLowEnergyLimit(0.1*keV);
    em->SetHighEnergyLimit(100.0*TeV);
    AddEmModel(1, em, fm);
  }
}
