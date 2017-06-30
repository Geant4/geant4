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
/// \file processes/phonon/src/G4PhononReflection.cc
/// \brief Implementation of the G4PhononReflection class
//
// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly 
// geometry dependent.Currently, phonons are killed when
// they reach a boundary. If the other side of the 
// boundary was Al, a hit is registered.
//  
// $Id: G4PhononReflection.cc 76799 2013-11-15 20:30:53Z mkelsey $
//
// 20131115  Throw exception if track's polarization state is invalid.

#include "G4PhononReflection.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GeometryTolerance.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"


G4PhononReflection::G4PhononReflection(const G4String& aName)
  : G4VPhononProcess(aName),
    kCarTolerance(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {;}

G4PhononReflection::~G4PhononReflection() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Always return DBL_MAX and Forced. This ensures that the process is
// called at the end of every step. In PostStepDoIt the process
// decides whether the step encountered a volume boundary and a
// reflection should be applied

G4double G4PhononReflection::GetMeanFreePath(const G4Track&, G4double,
					     G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This process handles the interaction of phonons with
// boundaries. Implementation of this class is highly geometry
// dependent.Currently, phonons are killed when they reach a
// boundary. If the other side of the boundary was Al, a hit is
// registered.
  
G4VParticleChange* G4PhononReflection::PostStepDoIt(const G4Track& aTrack,
						    const G4Step& aStep) { 
  aParticleChange.Initialize(aTrack);
   
  //Check if current step is limited by a volume boundary
  G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
  if (postStepPoint->GetStepStatus()!=fGeomBoundary) {
    //make sure that correct phonon velocity is used after the step
    int pol = GetPolarization(aTrack);
    if (pol < 0 || pol > 2) {
      G4Exception("G4PhononReflection::PostStepDoIt","Phonon001",
		  EventMustBeAborted, "Track is not a phonon");
      return &aParticleChange;		// NOTE: Will never get here
    }

    // FIXME:  This should be using wave-vector, shouldn't it?
    G4double vg = theLattice->MapKtoV(pol, aTrack.GetMomentumDirection());
    
    //Since step was not a volume boundary, just set correct phonon velocity and return
    aParticleChange.ProposeVelocity(vg);
    return &aParticleChange;
  }
  
  // do nothing but return is the step is too short
  // This is to allow actual reflection where after
  // the first boundary crossing a second, infinitesimal
  // step occurs crossing back into the original volume
  if (aTrack.GetStepLength()<=kCarTolerance/2) { 
    return &aParticleChange;
  }
  
  G4double eKin = aTrack.GetKineticEnergy();     
  aParticleChange.ProposeNonIonizingEnergyDeposit(eKin);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  
  return &aParticleChange; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


