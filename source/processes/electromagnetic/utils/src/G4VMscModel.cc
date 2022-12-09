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
// File name:     G4VMscModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.03.2008
//
// Modifications:
//
//
// Class Description:
//
// General interface to msc models

// -------------------------------------------------------------------
//

#include "G4VMscModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4TransportationManager.hh"
#include "G4LossTableManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::G4VMscModel(const G4String& nam):
  G4VEmModel(nam), 
  lambdalimit(1.*CLHEP::mm),
  geomMin(1.e-6*CLHEP::mm),
  geomMax(1.e50*CLHEP::mm),
  fDisplacement(0.,0.,0.),
  steppingAlgorithm(fUseSafety)
{
  dedx = 2.0*CLHEP::MeV*CLHEP::cm2/CLHEP::g;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::~G4VMscModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleChangeForMSC* 
G4VMscModel::GetParticleChangeForMSC(const G4ParticleDefinition* p)
{
  // recomputed for each new run
  if(nullptr == safetyHelper) {
    safetyHelper = G4TransportationManager::GetTransportationManager()
      ->GetSafetyHelper();
    safetyHelper->InitialiseHelper();
  }
  G4ParticleChangeForMSC* change = nullptr;
  if (nullptr != pParticleChange) {
    change = static_cast<G4ParticleChangeForMSC*>(pParticleChange);
  } else {
    change = new G4ParticleChangeForMSC();
  }
  if(IsMaster() && nullptr != p) {

    // table is always built for low mass particles 
    if(p->GetParticleName() != "GenericIon" &&
       (p->GetPDGMass() < CLHEP::GeV || ForceBuildTableFlag()) ) {

      G4EmParameters* param = G4EmParameters::Instance();
      G4LossTableBuilder* builder = 
	G4LossTableManager::Instance()->GetTableBuilder();
      G4double emin = std::max(LowEnergyLimit(), LowEnergyActivationLimit());
      G4double emax = std::min(HighEnergyLimit(), HighEnergyActivationLimit());
      emin = std::max(emin, param->MinKinEnergy());
      emax = std::min(emax, param->MaxKinEnergy());
      if(emin < emax) {
	xSectionTable = builder->BuildTableForModel(xSectionTable, this, p, 
						    emin, emax, useSpline);
      }
    }
  }
  return change;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::InitialiseParameters(const G4ParticleDefinition* part)
{
  if(IsLocked()) { return; }
  G4EmParameters* param = G4EmParameters::Instance();
  if(std::abs(part->GetPDGEncoding()) == 11) {
    steppingAlgorithm = param->MscStepLimitType(); 
    facrange = param->MscRangeFactor(); 
    latDisplasment = param->LateralDisplacement();
  } else {
    steppingAlgorithm = param->MscMuHadStepLimitType(); 
    facrange = param->MscMuHadRangeFactor(); 
    latDisplasment = param->MuHadLateralDisplacement();
  }
  skin = param->MscSkin();
  facgeom = param->MscGeomFactor();
  facsafety = param->MscSafetyFactor();
  lambdalimit = param->MscLambdaLimit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::DumpParameters(std::ostream& out) const
{
  G4String alg = "UseSafety";
  if (steppingAlgorithm == fUseDistanceToBoundary) alg = "DistanceToBoundary";
  else if (steppingAlgorithm == fMinimal) alg = "Minimal";
  else if (steppingAlgorithm == fUseSafetyPlus) alg = "SafetyPlus";

  out << std::setw(18) << "StepLim=" << alg << " Rfact=" << facrange 
      << " Gfact=" << facgeom << " Sfact=" << facsafety << " DispFlag:" << latDisplasment
      << " Skin=" << skin << " Llim=" << lambdalimit/CLHEP::mm << " mm" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                    const G4MaterialCutsCouple*,
                                    const G4DynamicParticle*,
                                    G4double, G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                  

G4double G4VMscModel::GetDEDX(const G4ParticleDefinition* part, G4double kinEnergy,
                              const G4MaterialCutsCouple* couple)
{
  G4double x;
  if (nullptr != ionisation) {
    x = ionisation->GetDEDX(kinEnergy, couple);
  } else {
    const G4double q = part->GetPDGCharge()*inveplus;
    x = dedx*q*q;
  }
  return x;
}                                                                                                                                                                                                                 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                  

G4double G4VMscModel::GetDEDX(const G4ParticleDefinition* part, G4double kinEnergy,
                              const G4MaterialCutsCouple* couple, G4double logKinEnergy)
{
  G4double x;
  if (nullptr != ionisation) {
    x = ionisation->GetDEDX(kinEnergy, couple, logKinEnergy);
  } else {
    const G4double q = part->GetPDGCharge()*inveplus;
    x = dedx*q*q;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMscModel::GetRange(const G4ParticleDefinition* part, G4double kinEnergy,
                               const G4MaterialCutsCouple* couple)
{
  //  << ionisation << "  " << part->GetParticleName() << G4endl;
  localtkin = kinEnergy;
  if (nullptr != ionisation) {
    localrange = ionisation->GetRange(kinEnergy, couple);
  } else {
    const G4double q = part->GetPDGCharge()*inveplus;
    localrange = kinEnergy/(dedx*q*q*couple->GetMaterial()->GetDensity());
  }
  //G4cout << "R(mm)= " << localrange << "  "  << ionisation << G4endl;
  return localrange;
}                                                                                                                                                                                                                 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                  

G4double G4VMscModel::GetRange(const G4ParticleDefinition* part,G4double kinEnergy,
                               const G4MaterialCutsCouple* couple, G4double logKinEnergy)
{
  //G4cout << "G4VMscModel::GetRange E(MeV)= " << kinEnergy << "  "
  //  << ionisation << "  " << part->GetParticleName() << G4endl;
  localtkin = kinEnergy;
  if (nullptr != ionisation) {
    localrange = ionisation->GetRange(kinEnergy, couple, logKinEnergy);
  } else {
    const G4double q = part->GetPDGCharge()*inveplus;
    localrange = kinEnergy/(dedx*q*q*couple->GetMaterial()->GetDensity());
  }
  //G4cout << "R(mm)= " << localrange << "  "  << ionisation << G4endl;
  return localrange;
}                                                                                                                                                                                                                 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMscModel::GetEnergy(const G4ParticleDefinition* part,
                                G4double range, const G4MaterialCutsCouple* couple)
{
  G4double e;
  //G4cout << "G4VMscModel::GetEnergy R(mm)= " << range << "  " << ionisation
  //     << "  Rlocal(mm)= " << localrange << "  Elocal(MeV)= " << localtkin << G4endl;
  if(nullptr != ionisation) { e = ionisation->GetKineticEnergy(range, couple); }
  else {
    e = localtkin;
    if(localrange > range) {
      G4double q = part->GetPDGCharge()*inveplus;
      e -= (localrange - range)*dedx*q*q*couple->GetMaterial()->GetDensity();
    }
  }
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
