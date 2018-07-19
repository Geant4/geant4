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
/// \file processes/phonon/src/G4PhononScattering.cc
/// \brief Implementation of the G4PhononScattering class
//
// $Id: G4PhononScattering.cc 76499 2013-11-12 05:33:22Z mkelsey $
//
// 20131111  Add verbose output for MFP calculation

#include "G4PhononScattering.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononLong.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"


G4PhononScattering::G4PhononScattering(const G4String& aName)
  : G4VPhononProcess(aName) {;}

G4PhononScattering::~G4PhononScattering() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PhononScattering::GetMeanFreePath(const G4Track& aTrack,
					     G4double /*previousStepSize*/,
					     G4ForceCondition* condition) {
  //Dynamical constants retrieved from PhysicalLattice
  G4double B = theLattice->GetScatteringConstant();
  G4double Eoverh = aTrack.GetKineticEnergy()/h_Planck;

  //Calculate mean free path
  G4double mfp = aTrack.GetVelocity()/(Eoverh*Eoverh*Eoverh*Eoverh*B);

  if (verboseLevel > 1)
    G4cout << "G4PhononScattering::GetMeanFreePath = " << mfp << G4endl;

  *condition = NotForced;
 
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PhononScattering::PostStepDoIt( const G4Track& aTrack,
						     const G4Step& aStep) {
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()==fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }
  
  //Initialize particle change
  aParticleChange.Initialize(aTrack);
  
  //randomly generate a new direction and polarization state
  G4ThreeVector newDir = G4RandomDirection();
  G4int polarization = ChoosePolarization(theLattice->GetLDOS(),
					  theLattice->GetSTDOS(),
					  theLattice->GetFTDOS());

  // Generate the new track after scattering
  // FIXME:  If polarization state is the same, just step the track!
  G4Track* sec =
    CreateSecondary(polarization, newDir, aTrack.GetKineticEnergy());
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(sec);

  // Scattered phonon replaces current track
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
