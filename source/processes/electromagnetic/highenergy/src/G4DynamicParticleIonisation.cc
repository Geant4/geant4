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
// File name:     G4DynamicParticleIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17.08.2024
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DynamicParticleIonisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticleFluctuation.hh"
#include "G4EmSecondaryParticleType.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4EmProcessSubType.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Log.hh"

namespace
{
  constexpr G4double ekinLimit = 0.2*CLHEP::MeV;
  const G4double twoln10 = 2*G4Log(10.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticleIonisation::G4DynamicParticleIonisation()
  : G4VContinuousDiscreteProcess("dynPartIoni")
{
  SetVerboseLevel(1);
  SetProcessSubType(fDynamicIonisation);
  theElectron = G4Electron::Electron();
 
  lManager = G4LossTableManager::Instance();  
  lManager->Register(this);

  fUrban = new G4DynamicParticleFluctuation();
  
  // define these flags only once
  auto param = G4EmParameters::Instance();
  fFluct = param->LossFluctuation();
  fLinLimit = 5*param->LinearLossLimit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticleIonisation::~G4DynamicParticleIonisation()
{
  lManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DynamicParticleIonisation::BuildPhysicsTable(const G4ParticleDefinition&)
{
  auto theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable(); 
  fCuts = theCoupleTable->GetEnergyCutsVector(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DynamicParticleIonisation::PreStepInitialisation(const G4Track& track)
{
  fCouple = track.GetMaterialCutsCouple();
  fMaterial = fCouple->GetMaterial();
  auto dpart = track.GetDynamicParticle();
  fEkinPreStep = dpart->GetKineticEnergy();
  fMass = std::max(dpart->GetMass(), CLHEP::electron_mass_c2);
  fCharge = dpart->GetCharge()/CLHEP::eplus;
  fRatio = fMass/CLHEP::proton_mass_c2;
  fLowestEkin = ekinLimit*fRatio;
  G4double tau = fEkinPreStep/fMass;
  G4double ratio = CLHEP::electron_mass_c2/fMass;
  fTmax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.) /
    (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  fCut = (*fCuts)[fCouple->GetIndex()];
  fCut = std::max(fCut, fMaterial->GetIonisation()->GetMeanExcitationEnergy());
  fCut = std::min(fCut, fTmax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleIonisation::AlongStepGetPhysicalInteractionLength(
                                  const G4Track&, G4double, G4double, G4double&,
                                  G4GPILSelection* selection)
{
  *selection = CandidateForSelection;

  // no step limit for the time being
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleIonisation::PostStepGetPhysicalInteractionLength(
                                  const G4Track& track, G4double previousStepSize,
                                  G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double x = DBL_MAX;
  G4double xsec = 0.0;

  PreStepInitialisation(track);

  if (fCharge != 0.0) {
    xsec = ComputeCrossSection(fEkinPreStep);
  }

  if (xsec <= 0.0) {
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;

  } else {
    if (theNumberOfInteractionLengthLeft < 0.0) {

      theNumberOfInteractionLengthLeft = -G4Log( G4UniformRand() );
      theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 

    } else if(currentInteractionLength < DBL_MAX) {

      // subtract NumberOfInteractionLengthLeft using previous step
      theNumberOfInteractionLengthLeft -= 
        previousStepSize/currentInteractionLength;

      theNumberOfInteractionLengthLeft = 
        std::max(theNumberOfInteractionLengthLeft, 0.0);
    }
    currentInteractionLength = 1.0/xsec;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;
  }
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4DynamicParticleIonisation::PostStepGetPhysicalInteractionLength ";
    G4cout << "  Process: " << GetProcessName()
           << " for unknown particle Mass(GeV)=" << fMass/CLHEP::GeV 
           << " charge=" << fCharge 
           << "  Material " << fMaterial->GetName()
           << "  Ekin(MeV)=" << fEkinPreStep/CLHEP::MeV 
           << "  MFP(cm)=" << currentInteractionLength/CLHEP::cm 
           << "  ProposedLength(cm)=" << x/CLHEP::cm <<G4endl;
  }
#endif
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
G4DynamicParticleIonisation::AlongStepDoIt(const G4Track& track,
					   const G4Step& step)
{
  fParticleChange.InitializeForAlongStep(track);

  // no energy loss
  if (fCharge == 0.0) { return &fParticleChange; }

  // stop low-energy object
  if (fEkinPreStep <= fLowestEkin) {
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(fEkinPreStep);
    return &fParticleChange;
  }
  
  G4double length = step.GetStepLength();
  G4double dedxPre = ComputeDEDX(fEkinPreStep);
  G4double eloss = dedxPre*length;
  G4double ekinPostStep = fEkinPreStep - eloss;

  // correction for large step if it is not the last step
  if (fEkinPreStep*fLinLimit < eloss && ekinPostStep > fLowestEkin) {
    G4double dedxPost = ComputeDEDX(ekinPostStep);
    eloss = (eloss + dedxPost*length)*0.5;
  }

  // do not sample fluctuations at the last step
  if (fFluct && fEkinPreStep > eloss) {
    eloss = fUrban->SampleFluctuations(fCouple, track.GetDynamicParticle(),
				       fCut, fTmax, length, eloss);
  }

  ekinPostStep = fEkinPreStep - eloss;

  // stop low-energy object
  if (ekinPostStep <= fLowestEkin) {
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(fEkinPreStep);
  } else {
    fParticleChange.SetProposedKineticEnergy(ekinPostStep);
    fParticleChange.ProposeLocalEnergyDeposit(eloss);
  }
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
G4DynamicParticleIonisation::PostStepDoIt(const G4Track& track, const G4Step&)
{
  theNumberOfInteractionLengthLeft = -1.0;
  fParticleChange.InitializeForPostStep(track);

  auto dp = track.GetDynamicParticle();
  G4double kinEnergy = dp->GetKineticEnergy();
  const G4double totEnergy = kinEnergy + fMass;
  const G4double beta2 = kinEnergy*(kinEnergy + 2.0*fMass)/(totEnergy*totEnergy);

  G4double deltaKinEnergy, f; 

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4double rndm[2];

  // sampling without nuclear size effect
  do {
    rndmEngineMod->flatArray(2, rndm);
    deltaKinEnergy = fCut*fTmax/(fCut*(1.0 - rndm[0]) + fTmax*rndm[0]);
    f = 1.0 - beta2*deltaKinEnergy/fTmax;
    // Loop checking, 14-Aug-2024, Vladimir Ivanchenko
  } while( rndm[1] > f);

  G4double deltaMomentum =
    std::sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*CLHEP::electron_mass_c2));
  G4double cost = deltaKinEnergy * (totEnergy + CLHEP::electron_mass_c2) /
      (deltaMomentum * dp->GetTotalMomentum());
  cost = std::min(cost, 1.0);
  const G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  const G4double phi = CLHEP::twopi*rndmEngineMod->flat();

  G4ThreeVector deltaDirection(sint*std::cos(phi), sint*std::sin(phi), cost);
  deltaDirection.rotateUz(dp->GetMomentumDirection());
  
  // create G4DynamicParticle object for delta ray
  auto delta = new G4DynamicParticle(theElectron, deltaDirection, deltaKinEnergy);
  auto t = new G4Track(delta, track.GetGlobalTime(), track.GetPosition());
  t->SetTouchableHandle(track.GetTouchableHandle());
  t->SetCreatorModelID(fSecID);
  fParticleChange.AddSecondary(t);

  // Change kinematics of primary particle
  kinEnergy -= deltaKinEnergy;
  G4ThreeVector finalP = dp->GetMomentum() - delta->GetMomentum();
  finalP = finalP.unit();
  
  fParticleChange.SetProposedKineticEnergy(kinEnergy);
  fParticleChange.SetProposedMomentumDirection(finalP);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleIonisation::ComputeDEDX(G4double ekin)
{
  G4double tau = ekin/fMass;
  G4double gam = tau + 1.0;
  G4double bg2 = tau * (tau + 2.0);
  G4double beta2 = bg2/(gam*gam);
  G4double xc = fCut/fTmax;

  G4double exc  = fMaterial->GetIonisation()->GetMeanExcitationEnergy();
  G4double exc2 = exc*exc;

  // general Bethe-Bloch formula
  G4double dedx = G4Log(2.0*CLHEP::electron_mass_c2*bg2*fCut/exc2) - (1.0 + xc)*beta2;

  // density correction
  G4double x = G4Log(bg2)/twoln10;
  dedx -= fMaterial->GetIonisation()->DensityCorrection(x);

  // now compute the total ionization loss per volume
  dedx *= CLHEP::twopi_mc2_rcl2*fCharge*fCharge*fMaterial->GetElectronDensity()/beta2;
  dedx = std::max(dedx, 0.0);
  return dedx; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleIonisation::ComputeCrossSection(G4double ekin)
{
  G4double cross = 0.0;
  if (fCut < fTmax) {

    G4double totEnergy = ekin + fMass;
    G4double energy2 = totEnergy*totEnergy;
    G4double beta2 = ekin*(ekin + 2.0*fMass)/energy2;

    cross = (fTmax - fCut)/(fCut*fTmax*beta2) - G4Log(fTmax/fCut)/fTmax;
    cross *= CLHEP::twopi_mc2_rcl2*fCharge*fCharge*fMaterial->GetElectronDensity();
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleIonisation::GetMeanFreePath(const G4Track& /* track */, G4double,
                                                      G4ForceCondition* condition)
{
  // Note: this method is not used at run-time, so its implementation is simplified.
  //       It might be eventually refined later.
  *condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DynamicParticleIonisation::GetContinuousStepLimit(const G4Track&, G4double,
							     G4double, G4double&)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DynamicParticleIonisation::ProcessDescription(std::ostream& out) const
{
  out << "G4DynamicParticleIonisation: dynamic ionisation" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
