// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleChangeForMSC.cc,v 1.1 1999-01-07 16:14:26 gunter Exp $
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
    G4cerr << "G4ParticleChangeForMSC::G4ParticleChangeForMSC() " << endl;
  }
#endif
}

G4ParticleChangeForMSC::~G4ParticleChangeForMSC() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cerr << "G4ParticleChangeForMSC::~G4ParticleChangeForMSC() " << endl;
  }
#endif
}

G4ParticleChangeForMSC::G4ParticleChangeForMSC(
             const G4ParticleChangeForMSC &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cerr << "G4ParticleChangeForMSC::  copy constructor is called " << endl;
   }
   *this = right;
}

// assignment operator
G4ParticleChangeForMSC & G4ParticleChangeForMSC::operator=(
                                   const G4ParticleChangeForMSC &right)
{
   if (verboseLevel>1) {
    G4cerr << "G4ParticleChangeForMSC:: assignment operator is called " << endl;
   }
   if (this != &right)
   {
      theMomentumChange = right.theMomentumChange;
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
  theMomentumChange        = pParticle->GetMomentumDirection();

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
  pPostStepPoint->SetMomentumDirection(theMomentumChange);

  // update position 
  pPostStepPoint->SetPosition( thePositionChange  );

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
       << setw(20) << thePositionChange.x()/mm
       << endl; 
  G4cout << "        Position - y (mm)   : " 
       << setw(20) << thePositionChange.y()/mm
       << endl; 
  G4cout << "        Position - z (mm)   : " 
       << setw(20) << thePositionChange.z()/mm
       << endl;
  G4cout << "        Momentum Direct - x : " 
       << setw(20) << theMomentumChange.x()
       << endl;
  G4cout << "        Momentum Direct - y : " 
       << setw(20) << theMomentumChange.y()
       << endl;
  G4cout << "        Momentum Direct - z : " 
       << setw(20) << theMomentumChange.z()
       << endl;
}







