// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleChangeForMSC.cc,v 1.6 2000-02-13 15:08:28 kurasige Exp $
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

#include "G4ParticleChangeForMSC.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

G4ParticleChangeForMSC::G4ParticleChangeForMSC():G4VParticleChange()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForMSC::G4ParticleChangeForMSC() " << G4endl;
  }
#endif
}

G4ParticleChangeForMSC::~G4ParticleChangeForMSC() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForMSC::~G4ParticleChangeForMSC() " << G4endl;
  }
#endif
}

G4ParticleChangeForMSC::G4ParticleChangeForMSC(
             const G4ParticleChangeForMSC &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForMSC::  copy constructor is called " << G4endl;
   }
      theMomentumDirectionChange = right.theMomentumDirectionChange;
      thePositionChange = right.thePositionChange;
      theTrueStepLength = right.theTrueStepLength;
}

// assignment operator
G4ParticleChangeForMSC & G4ParticleChangeForMSC::operator=(
                                   const G4ParticleChangeForMSC &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForMSC:: assignment operator is called " << G4endl;
   }
   if (this != &right)
   {  
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;


      theMomentumDirectionChange = right.theMomentumDirectionChange;
      thePositionChange = right.thePositionChange;
      theTrueStepLength = right.theTrueStepLength;
   }
   return *this;
}

//----------------------------------------------------------------
// functions for Initialization
//

void G4ParticleChangeForMSC::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  // set Energy/Momentum etc. equal to those of the parent particle
  const G4DynamicParticle*  pParticle = track.GetDynamicParticle();
  theMomentumDirectionChange        = pParticle->GetMomentumDirection();

  // set Position equal to those of the parent track
  thePositionChange      = track.GetPosition();
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4ParticleChangeForMSC::UpdateStepForAlongStep(G4Step* pStep)
{
  
  //  Update the G4Step specific attributes 
  pStep->SetStepLength(theTrueStepLength) ;

  return pStep;
}

G4Step* G4ParticleChangeForMSC::UpdateStepForPostStep(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  // Take note that the return type of GetMomentumChange is a
  // pointer to G4ParticleMometum. Also it is a normalized 
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
 
  // update  momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);

  // update position 
  pPostStepPoint->SetPosition( thePositionChange  );

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

G4Step* G4ParticleChangeForMSC::UpdateStepForAtRest(G4Step* pStep)
{ 

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

//----------------------------------------------------------------
// methods for printing messages  
//

void G4ParticleChangeForMSC::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4cout.precision(3);
  G4cout << "        Position - x (mm)   : " 
       << G4std::setw(20) << thePositionChange.x()/mm
       << G4endl; 
  G4cout << "        Position - y (mm)   : " 
       << G4std::setw(20) << thePositionChange.y()/mm
       << G4endl; 
  G4cout << "        Position - z (mm)   : " 
       << G4std::setw(20) << thePositionChange.z()/mm
       << G4endl;
  G4cout << "     Momentum Direction - x : " 
       << G4std::setw(20) << theMomentumDirectionChange.x()
       << G4endl;
  G4cout << "     Momentum Direction - y : " 
       << G4std::setw(20) << theMomentumDirectionChange.y()
       << G4endl;
  G4cout << "     Momentum Direction - z : " 
       << G4std::setw(20) << theMomentumDirectionChange.z()
       << G4endl;
}


G4bool G4ParticleChangeForMSC::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // check     

  // MomentumDirection should be unit vector
  accuracy = abs(theMomentumDirectionChange.mag2()-1.0);
  if (accuracy > accuracyForWarning) {
    G4cout << "  G4ParticleChangeForMSC::CheckIt  : ";
    G4cout << "the Momentum Change is not unit vector !!" << G4endl;
    G4cout << "  Difference:  " << accuracy << G4endl;
    itsOK = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // dump out information of this particle change
  if (!itsOK) { 
    G4cout << " G4ParticleChangeForMSC::CheckIt " <<G4endl;
    DumpInfo();
  }

  // Exit with error
  if (exitWithError) G4Exception("G4ParticleChangeForMSC::CheckIt");

  //correction
  if (!itsOK) {
    G4double vmag = theMomentumDirectionChange.mag();
    theMomentumDirectionChange = (1./vmag)*theMomentumDirectionChange;
  }


  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}




