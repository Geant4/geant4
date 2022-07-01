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
// File name:     G4MuPairProduction
//
// Author:        Laszlo Urban
//
// Creation date: 30.05.1998
//
// Modifications:
//
// 04-06-98 in DoIt,secondary production condition:
//          range>std::min(threshold,safety)
// 26-10-98 new stuff from R. Kokoulin + cleanup , L.Urban
// 06-05-99 bug fixed , L.Urban
// 10-02-00 modifications+bug fix , new e.m. structure, L.Urban
// 29-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 20-09-01 (L.Urban) in ComputeMicroscopicCrossSection, remove:
//          if(MaxPairEnergy<CutInPairEnergy) MaxPairEnergy=CutInPairEnergy
// 26-09-01 completion of store/retrieve PhysicsTable
// 28-09-01 suppression of theMuonPlus ..etc..data members (mma)
// 29-10-01 all static functions no more inlined (mma)
// 07-11-01 particleMass becomes a local variable (mma)
// 19-08-02 V.Ivanchenko update to new design
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 27-09-03 e+ set to be a secondary particle (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 10-02-04 Add lowestKinEnergy (V.Ivanchenko)
// 17-08-04 Utilise mu+ tables for mu- (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuPairProduction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"
#include "G4VEmModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4ElementData.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuPairProduction::G4MuPairProduction(const G4String& name)
  : G4VEnergyLossProcess(name),
    lowestKinEnergy(0.85*CLHEP::GeV)
{
  SetProcessSubType(fPairProdByCharged);
  SetSecondaryParticle(G4Positron::Positron());
  SetIonisation(false);
  SetSpline(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuPairProduction::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProduction::MinPrimaryEnergy(const G4ParticleDefinition*,
					      const G4Material*,
					      G4double)
{
  return lowestKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuPairProduction::InitialiseEnergyLossProcess(
                         const G4ParticleDefinition* part,
			 const G4ParticleDefinition*)
{
  if (isInitialised) { return; }
  isInitialised = true;

  theParticle = part;

  G4VEmModel* mod = EmModel(0);
  if(nullptr == mod) {
    lowestKinEnergy = std::max(lowestKinEnergy, part->GetPDGMass()*8.0);
    auto ptr = new G4MuPairProductionModel(part);
    ptr->SetLowestKineticEnergy(lowestKinEnergy);
    mod = ptr; 
    SetEmModel(mod);
  }

  G4VEmFluctuationModel* fm = nullptr;
  G4EmParameters* param = G4EmParameters::Instance();
  mod->SetLowEnergyLimit(param->MinKinEnergy());
  mod->SetHighEnergyLimit(param->MaxKinEnergy());
  mod->SetSecondaryThreshold(param->MuHadBremsstrahlungTh());
  AddEmModel(1, mod, fm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuPairProduction::StreamProcessInfo(std::ostream& out) const
{
  G4ElementData* ed = EmModel()->GetElementData();
  if(ed) {
    for(G4int Z=1; Z<93; ++Z) {
      G4Physics2DVector* pv = ed->GetElement2DData(Z);
      if(pv) {
        out << "      Sampling table " << pv->GetLengthY()
	    << "x" << pv->GetLengthX() << "; from "
	    << std::exp(pv->GetY(0))/GeV << " GeV to " 
	    << std::exp(pv->GetY(pv->GetLengthY()-1))/TeV 
	    << " TeV " << G4endl;
	break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuPairProduction::ProcessDescription(std::ostream& out) const
{
  out << "  Electron-positron pair production by muons";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
