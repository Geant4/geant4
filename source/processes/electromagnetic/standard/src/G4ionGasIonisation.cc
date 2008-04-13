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
// $Id: G4ionGasIonisation.cc,v 1.8 2008-04-13 18:06:37 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ionGasIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.07.2007
//
// Modifications:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionGasIonisation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4ionGasIonisation::G4ionGasIonisation(const G4String& name)
  : G4ionIonisation(name),
    currParticle(0),
    baseParticle(0),
    initialised(false)
{
  atomXS = CLHEP::pi*CLHEP::Bohr_radius*CLHEP::Bohr_radius;
  verboseLevel = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionGasIonisation::~G4ionGasIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionGasIonisation::InitialiseEnergyLossProcess(
		      const G4ParticleDefinition* part,
		      const G4ParticleDefinition* bpart)
{
  G4ionIonisation::InitialiseEnergyLossProcess(part, bpart);
  currTh = BetheBlochEnergyThreshold();
  if(initialised) return;

  currParticle = part;

  if(part == bpart || part == G4GenericIon::GenericIon()) baseParticle = 0;
  else if(bpart == 0) baseParticle = G4GenericIon::GenericIon();
  else                baseParticle = bpart;

  if(baseParticle) basePartMass = baseParticle->GetPDGMass();
  else             basePartMass = currParticle->GetPDGMass();
  
  initialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionGasIonisation::PrintInfo()
{
  G4ionIonisation::PrintInfo();
  G4cout << "      Version of ion process with simulation discrete ion/media change exchange."
	 << G4endl;	 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionGasIonisation::InitialiseMassCharge(const G4Track& track)
{
  // First step of an ion
  if(track.GetCurrentStepNumber() == 1) {
    currParticle = track.GetDefinition();
    ionZ = G4int(currParticle->GetPDGCharge()/eplus + 0.5);
    currentIonZ = G4int(track.GetDynamicParticle()->GetCharge()/eplus + 0.5);
    currMassRatio = basePartMass/currParticle->GetPDGMass();
  }
  // any step 
  G4double q = eplus*currentIonZ;
  currCharge2= currentIonZ*currentIonZ;
  SetDynamicMassCharge(currMassRatio, q*q);
  preStepKinEnergy = track.GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionGasIonisation::CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* dp,
					      G4double& eloss,
					      G4double& niel,
                                              G4double s)
{
  const G4ParticleDefinition* part = dp->GetDefinition();
  const G4Material* mat = couple->GetMaterial();

  // estimate mean energy at the step
  G4double e = preStepKinEnergy - eloss*0.5;
  if(e <= 0.0) e = 0.5*preStepKinEnergy;

  // add corrections
  if(eloss < preStepKinEnergy) {

    // use Bethe-Bloch with corrections
    if(e*currMassRatio > currTh) {
      eloss += s*corr->IonHighOrderCorrections(part,mat,e);

      // Correction for data points 
    } else {
      //      if(stopDataActive)
      eloss *= corr->EffectiveChargeCorrection(part,mat,e)*currCharge2
	/effCharge->EffectiveChargeSquareRatio(part,mat,e);
      // G4cout<<"Below th: eloss= "<<eloss<<" f= "<<eloss/eloss0<<G4endl;
    }

    // effective number of collisions
    G4double x = mat->GetElectronDensity()*s*atomXS;
    // equilibrium charge
    G4double q = fParticleChange.GetProposedCharge(); 
  
    // sample charge change during the step
    fParticleChange.SetProposedCharge(SampleChargeAfterStep(q, x));
  }

  // use nuclear stopping 
  if(NuclearStoppingFlag() && 
     preStepKinEnergy*currMassRatio < 50.*currTh*currCharge2) {

    G4double nloss = s*corr->NuclearDEDX(part,mat,e);
    if(eloss + nloss > preStepKinEnergy) {
      nloss *= (preStepKinEnergy/(eloss + nloss));
      eloss = preStepKinEnergy - nloss;
    }
    niel += nloss;
    //    G4cout << "G4ionIonisation::CorrectionsAlongStep: e= " << preKinEnergy
    //	   << " de= " << eloss << " NIEL= " << nloss << G4endl;
    fParticleChange.ProposeNonIonizingEnergyDeposit(nloss);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionGasIonisation::SampleChargeAfterStep(G4double qeff, G4double xeff)
{
  // qeff - equilibrium charge
  // xeff - effective number of collisions
  // q    - current charge
  G4double q = eplus*currentIonZ;
  if(verboseLevel > 1) G4cout << "G4ionGasIonisation: Q1= " << currentIonZ
			      << " Qeff= " << qeff/eplus << "  Neff= " << xeff
			      << G4endl;
  return q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
