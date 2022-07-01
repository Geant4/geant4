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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eBremsstrahlung
//
// Author:        Michel Maire
//
// Creation date: 26.06.1996
//
// Modified by Michel Maire, Vladimir Ivanchenko, Andreas Schaelicke  
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eBremsstrahlung.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;
 
G4eBremsstrahlung::G4eBremsstrahlung(const G4String& name):
  G4VEnergyLossProcess(name)
{
  SetProcessSubType(fBremsstrahlung);
  SetSecondaryParticle(G4Gamma::Gamma());
  SetIonisation(false);
  SetCrossSectionType(fEmTwoPeaks);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlung::~G4eBremsstrahlung() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4eBremsstrahlung::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4eBremsstrahlung::InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					       const G4ParticleDefinition*)
{
  if(!isInitialised) {
    G4EmParameters* param = G4EmParameters::Instance();

    G4double emax = param->MaxKinEnergy();
    G4VEmFluctuationModel* fm = nullptr;

    if (nullptr == EmModel(0)) { SetEmModel(new G4SeltzerBergerModel()); }
    G4double energyLimit = std::min(EmModel(0)->HighEnergyLimit(), CLHEP::GeV);
    EmModel(0)->SetHighEnergyLimit(energyLimit);
    EmModel(0)->SetSecondaryThreshold(param->BremsstrahlungTh());
    EmModel(0)->SetLPMFlag(false);
    AddEmModel(1, EmModel(0), fm);

    if(emax > energyLimit) {
      if (nullptr == EmModel(1)) { 
	SetEmModel(new G4eBremsstrahlungRelModel());
      }
      EmModel(1)->SetLowEnergyLimit(energyLimit);
      EmModel(1)->SetHighEnergyLimit(emax); 
      EmModel(1)->SetSecondaryThreshold(param->BremsstrahlungTh());
      EmModel(1)->SetLPMFlag(param->LPM());
      AddEmModel(1, EmModel(1), fm);
    }
    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlung::StreamProcessInfo(std::ostream& out) const
{
  if(nullptr != EmModel(0)) {
    G4EmParameters* param = G4EmParameters::Instance();
    G4double eth = param->BremsstrahlungTh(); 
    out << "      LPM flag: " << param->LPM() << " for E > " 
	<< EmModel(0)->HighEnergyLimit()/GeV << " GeV";
    if(eth < DBL_MAX) { 
      out << ",  VertexHighEnergyTh(GeV)= " << eth/GeV; 
    }
    out << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

void G4eBremsstrahlung::ProcessDescription(std::ostream& out) const
{
  out << "  Bremsstrahlung";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
