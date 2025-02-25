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
// File name:     G4DynamicParticleMSC
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17.08.2024
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DynamicParticleMSC.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmProcessSubType.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

namespace
{
  constexpr G4double c_highland = 13.6*CLHEP::MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticleMSC::G4DynamicParticleMSC()
  : G4VContinuousDiscreteProcess("dynPartMSC")
{
  SetVerboseLevel(1);
  SetProcessSubType(fDynamicMultipleScattering);
  lManager = G4LossTableManager::Instance();  
  lManager->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticleMSC::~G4DynamicParticleMSC()
{
  lManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DynamicParticleMSC::PreStepInitialisation(const G4Track& track)
{
  fMaterial = track.GetMaterial();
  fZeff = fMaterial->GetIonisation()->GetZeffective();
  auto dpart = track.GetDynamicParticle();
  fEkinPreStep = dpart->GetKineticEnergy();
  fBeta = dpart->GetBeta();
  fCharge = dpart->GetCharge()/CLHEP::eplus;
  fMass = std::max(dpart->GetMass(), CLHEP::electron_mass_c2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleMSC::AlongStepGetPhysicalInteractionLength(
                                  const G4Track& track, G4double, G4double, G4double&,
                                  G4GPILSelection* selection)
{
  *selection = CandidateForSelection;
  PreStepInitialisation(track);

  // no step limit for the time being
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleMSC::PostStepGetPhysicalInteractionLength(
                                  const G4Track&, G4double,
                                  G4ForceCondition* condition)
{
  *condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4DynamicParticleMSC::AlongStepDoIt(const G4Track& track,
						              const G4Step& step)
{
  fParticleChange.InitialiseMSC(track, step);

  // no energy loss
  if (fCharge == 0.0) { return &fParticleChange; }

  G4double geomLength = step.GetStepLength();
  G4double y = geomLength/fMaterial->GetRadlen();
  G4double theta0 = c_highland*std::abs(fCharge)*std::sqrt(y)*
    (1.0 + 0.038*G4Log(y*fCharge*fCharge/(fBeta*fBeta)))/fBeta;

  if (theta0 < 0.001) { return &fParticleChange; }
  G4double cost = 1.0;
  G4double r = G4UniformRand();
  if (theta0 < 1.0) {
    G4double theta2 = theta0*theta0;
    cost -= theta2*G4Log(1.0 + r*(G4Exp(2.0/theta2) - 1.0));
  } else {
    cost -= 2.0*r;
  }
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  fNewDir.set(sint*std::cos(phi), sint*std::sin(phi), cost);
  fNewDir.rotateUz(step.GetPostStepPoint()->GetMomentumDirection());
  
  fParticleChange.ProposeMomentumDirection(fNewDir);
  fParticleChange.ProposeTrueStepLength(geomLength);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleMSC::GetMeanFreePath(const G4Track&, G4double,
					       G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleMSC::GetContinuousStepLimit(const G4Track& track,
						      G4double previousStepSize,
                                                      G4double currentMinimalStep,
						      G4double& currentSafety)
{
  G4GPILSelection selection = NotCandidateForSelection;
  G4double x = AlongStepGetPhysicalInteractionLength(track, previousStepSize,
						     currentMinimalStep,
					             currentSafety, &selection);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DynamicParticleMSC::ProcessDescription(std::ostream& out) const
{
  out << "G4DynamicParticleMSC: no delta rays" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
