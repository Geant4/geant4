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
// $Id: G4VMscModel.cc 105745 2017-08-16 13:14:37Z gcosmo $
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
  safetyHelper(nullptr),
  ionisation(nullptr),
  facrange(0.04),
  facgeom(2.5),
  facsafety(0.3),
  skin(1.0),
  dtrl(0.05),
  lambdalimit(mm),
  geomMin(1.e-6*CLHEP::mm),
  geomMax(1.e50*CLHEP::mm),
  fDisplacement(0.,0.,0.),
  steppingAlgorithm(fUseSafety),
  samplez(false),
  latDisplasment(true)
{
  dedx       = 2.0*CLHEP::MeV*CLHEP::cm2/CLHEP::g;
  localrange = DBL_MAX;
  localtkin  = 0.0;
  currentPart = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::~G4VMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleChangeForMSC* 
G4VMscModel::GetParticleChangeForMSC(const G4ParticleDefinition* p)
{
  // recomputed for each new run
  if(!safetyHelper) {
    safetyHelper = G4TransportationManager::GetTransportationManager()
      ->GetSafetyHelper();
    safetyHelper->InitialiseHelper();
  }
  G4ParticleChangeForMSC* change = nullptr;
  if (pParticleChange) {
    change = static_cast<G4ParticleChangeForMSC*>(pParticleChange);
  } else {
    change = new G4ParticleChangeForMSC();
  }
  if(p) {

    // table is never built for GenericIon 
    if(p->GetParticleName() == "GenericIon") {
      if(xSectionTable) {
        xSectionTable->clearAndDestroy();
        delete xSectionTable;
        xSectionTable = nullptr;
      }

      // table is always built for low mass particles 
    } else if(p->GetPDGMass() < 4.5*GeV || ForceBuildTableFlag()) {

      G4EmParameters* param = G4EmParameters::Instance();
      idxTable = 0;
      G4LossTableBuilder* builder = 
	G4LossTableManager::Instance()->GetTableBuilder();
      if(IsMaster()) {
        G4double emin = std::max(LowEnergyLimit(), LowEnergyActivationLimit());
        G4double emax = std::min(HighEnergyLimit(), HighEnergyActivationLimit());
        emin = std::max(emin, param->MinKinEnergy());
        emax = std::min(emax, param->MaxKinEnergy());
        if(emin < emax) {
          xSectionTable = builder->BuildTableForModel(xSectionTable, this, p, 
                                                      emin, emax, true);
        }
      }
      theDensityFactor = builder->GetDensityFactors();
      theDensityIdx = builder->GetCoupleIndexes();
    }
  }
  return change;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector& 
G4VMscModel::SampleScattering(const G4ThreeVector&, G4double)
{
  return fDisplacement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMscModel::ComputeTruePathLengthLimit(const G4Track&, G4double&)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMscModel::ComputeGeomPathLength(G4double truePathLength)
{
  return truePathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMscModel::ComputeTrueStepLength(G4double geomPathLength)
{
  return geomPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                    const G4MaterialCutsCouple*,
                                    const G4DynamicParticle*,
                                    G4double, G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
