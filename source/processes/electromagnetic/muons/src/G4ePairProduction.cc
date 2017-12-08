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
// $Id: G4ePairProduction.cc 84953 2014-10-22 15:08:57Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ePairProduction
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17.03.2016
//
// Modifications:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ePairProduction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4ElementData.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4ePairProduction::G4ePairProduction(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(nullptr),
    lowestKinEnergy(100.*MeV),
    isInitialised(false)
{
  SetProcessSubType(fPairProdByCharged);
  SetSecondaryParticle(G4Positron::Positron());
  SetIonisation(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ePairProduction::~G4ePairProduction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ePairProduction::IsApplicable(const G4ParticleDefinition& p)
{
  return (G4Electron::Electron() == &p || G4Positron::Positron() == &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ePairProduction::MinPrimaryEnergy(const G4ParticleDefinition*,
					     const G4Material*,
					     G4double)
{
  return lowestKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePairProduction::InitialiseEnergyLossProcess(
                         const G4ParticleDefinition* part,
			 const G4ParticleDefinition*)
{
  if (!isInitialised) {
    isInitialised = true;

    theParticle = part;

    G4MuPairProductionModel* mod = new G4MuPairProductionModel(part); 
    SetEmModel(mod);

    lowestKinEnergy = std::max(lowestKinEnergy, part->GetPDGMass()*8.0);
    mod->SetLowestKineticEnergy(lowestKinEnergy);

    G4VEmFluctuationModel* fm = nullptr;
    G4EmParameters* param = G4EmParameters::Instance();
    mod->SetLowEnergyLimit(param->MinKinEnergy());
    mod->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, mod, fm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePairProduction::StreamProcessInfo(std::ostream& out,
                                          G4String endOfLine) const
{
  G4ElementData* ed = EmModel(0)->GetElementData();
  if(ed) {
    for(G4int Z=1; Z<93; ++Z) {
      G4Physics2DVector* pv = ed->GetElement2DData(Z);
      if(pv) {
        out << "      Sampling table " << pv->GetLengthY()
	    << "x" << pv->GetLengthX() << "; from "
	    << G4Exp(pv->GetY(0))/GeV << " GeV to " 
	    << G4Exp(pv->GetY(pv->GetLengthY()-1))/TeV 
	    << " TeV " << endOfLine;
	break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePairProduction::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Pair production</strong>";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
