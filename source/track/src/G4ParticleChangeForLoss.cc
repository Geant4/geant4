// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleChangeForLoss.cc,v 1.4 1999-12-15 14:53:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
// --------------------------------------------------------------

#include "G4ParticleChangeForLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

G4ParticleChangeForLoss::G4ParticleChangeForLoss():G4VParticleChange()
{
  debugFlag = false;
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForLoss::G4ParticleChangeForLoss() " << G4endl;
  }
#endif
}

G4ParticleChangeForLoss::~G4ParticleChangeForLoss() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForLoss::~G4ParticleChangeForLoss() " << G4endl;
  }
#endif
}

// copy constructor
G4ParticleChangeForLoss::G4ParticleChangeForLoss(const G4ParticleChangeForLoss &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForLoss::  copy constructor is called " << G4endl;
   }
   theEnergyChange = right.theEnergyChange;
}

// assignemnt operator
G4ParticleChangeForLoss & G4ParticleChangeForLoss::operator=(const G4ParticleChangeForLoss &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForLoss:: assignment operator is called " << G4endl;
   }
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theTrueStepLength = right.theTrueStepLength;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
     
      theEnergyChange = right.theEnergyChange;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit ;
   }
   return *this;
}


//----------------------------------------------------------------
// functions for Initialization
//

void G4ParticleChangeForLoss::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  // set Energy equal to those of the parent particle
  const G4DynamicParticle*  pParticle = track.GetDynamicParticle();
  theEnergyChange          = pParticle->GetKineticEnergy();
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4ParticleChangeForLoss::UpdateStepForAlongStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint). 
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint. 
  

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
 
  // calculate new kinetic energy
  G4double energy = pPostStepPoint->GetKineticEnergy() 
                    + (theEnergyChange - pPreStepPoint->GetKineticEnergy()); 

  // update kinetic energy and momentum direction
  if (energy > 0.0) {
    pPostStepPoint->SetKineticEnergy( energy );
  } else {
    // stop case
    pPostStepPoint->SetKineticEnergy(0.0);
  }

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

//----------------------------------------------------------------
// methods for printing messages  
//

void G4ParticleChangeForLoss::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4cout.precision(3);
  G4cout << "        Kinetic Energy (MeV): " 
       << G4std::setw(20) << theEnergyChange/MeV
       << G4endl;
}

G4bool G4ParticleChangeForLoss::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // Energy should not be lager than initial value
  accuracy = ( theEnergyChange - aTrack.GetKineticEnergy())/MeV;
  if (accuracy > accuracyForWarning) {
    G4cout << "  G4ParticleChangeForLoss::CheckIt    : ";
    G4cout << "the energy becoes larger than the initial value !!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
    itsOK = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // dump out information of this particle change
  if (!itsOK) { 
    G4cout << " G4ParticleChangeForLoss::CheckIt " <<G4endl;
    G4cout << " pointer : " << this <<G4endl ;
    DumpInfo();
  }

  // Exit with error
  if (exitWithError) G4Exception("G4ParticleChangeForLoss::CheckIt");

  //correction
  if (!itsOK) {
    theEnergyChange = aTrack.GetKineticEnergy();
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}












