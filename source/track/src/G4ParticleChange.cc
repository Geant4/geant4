// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleChange.cc,v 1.7 1999-10-05 06:48:35 kurasige Exp $
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
//   Change default debug flag to false             10 May. 1998  H.Kurahige
//   Add Track weight                               12 Nov. 1998  H.Kurashige
//   Activate CheckIt method for VERBOSE mode       14 Dec. 1998 H.Kurashige
//   Modified CheckIt method for time                9 Feb. 1999 H.Kurashige
// --------------------------------------------------------------

#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

G4bool G4ParticleChange::fUseEBForAll = false;

G4ParticleChange::G4ParticleChange():G4VParticleChange(false)
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChange::G4ParticleChange() " << endl;
  }
#endif
}

G4ParticleChange::G4ParticleChange(G4bool useEB):G4VParticleChange(useEB)
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChange::G4ParticleChange() " << endl;
  }
#endif
}

G4ParticleChange::~G4ParticleChange() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChange::~G4ParticleChange() " << endl;
  }
#endif
}

// copy constructor
G4ParticleChange::G4ParticleChange(const G4ParticleChange &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChange::  copy constructor is called " << endl;
   }
   theMomentumDirectionChange = right.theMomentumDirectionChange;
   thePolarizationChange = right.thePolarizationChange;
   thePositionChange = right.thePositionChange;
   theTimeChange = right.theTimeChange;
   theEnergyChange = right.theEnergyChange;
   theMassChange = right.theMassChange;
   theChargeChange = right.theChargeChange;
   theWeightChange = right.theWeightChange;
}

// assignemnt operator
G4ParticleChange & G4ParticleChange::operator=(const G4ParticleChange &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChange:: assignment operator is called " << endl;
   }
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theMomentumDirectionChange = right.theMomentumDirectionChange;
      thePolarizationChange = right.thePolarizationChange;
      thePositionChange = right.thePositionChange;
      theTimeChange = right.theTimeChange;
      theEnergyChange = right.theEnergyChange;
      theMassChange = right.theMassChange;
      theChargeChange = right.theChargeChange;
      theWeightChange = right.theWeightChange;
      theTrueStepLength = right.theTrueStepLength;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
   }
   return *this;
}

G4bool G4ParticleChange::operator==(const G4ParticleChange &right) const
{
   return ((G4VParticleChange *)this == (G4VParticleChange *) &right);
}

G4bool G4ParticleChange::operator!=(const G4ParticleChange &right) const
{
   return ((G4VParticleChange *)this != (G4VParticleChange *) &right);
}


//----------------------------------------------------------------
// methods for handling secondaries 
//

void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle, 
				    G4bool   IsGoodForTracking    )
{
  //  create track
  G4Track* aTrack = new G4Track(aParticle, theTimeChange, thePositionChange);

  // set IsGoodGorTrackingFlag
  if (IsGoodForTracking) aTrack->SetGoodForTrackingFlag();

  //   Touchable is a temporary object, so you cannot keep the pointer
  aTrack->SetTouchable(0);

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle, 
				    G4ThreeVector      newPosition,
				    G4bool   IsGoodForTracking    )
{
  //  create track
  G4Track*  aTrack = new G4Track(aParticle, theTimeChange, newPosition);

  // set IsGoodGorTrackingFlag
  if (IsGoodForTracking) aTrack->SetGoodForTrackingFlag();

  //   Touchable is a temporary object, so you cannot keep the pointer
  aTrack->SetTouchable(0);

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle, 
				    G4double           newTime,
				    G4bool   IsGoodForTracking    )
{
  //  create track
  G4Track*  aTrack = new G4Track(aParticle, newTime, thePositionChange); 

  // set IsGoodGorTrackingFlag
  if (IsGoodForTracking) aTrack->SetGoodForTrackingFlag();
 
  //   Touchable is a temporary object, so you cannot keep the pointer
  aTrack->SetTouchable(0);

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChange::AddSecondary(G4Track* aTrack)
{
  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

//----------------------------------------------------------------
// functions for Initialization
//

void G4ParticleChange::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  // set Energy/Momentum etc. equal to those of the parent particle
  const G4DynamicParticle*  pParticle = track.GetDynamicParticle();
  theEnergyChange          = pParticle->GetKineticEnergy();
  theMomentumDirectionChange        = pParticle->GetMomentumDirection();
  thePolarizationChange    = pParticle->GetPolarization();
  theProperTimeChange      = pParticle->GetProperTime();

  // Set mass/charge of DynamicParticle
  theMassChange = pParticle->GetMass();
  theChargeChange = pParticle->GetCharge();

  // set Position/Time etc. equal to those of the parent track
  thePositionChange      = track.GetPosition();
  theTimeChange          = track.GetGlobalTime();

  theWeightChange        = track.GetWeight();
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4ParticleChange::UpdateStepForAlongStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint). 
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint. 
  
  // Take note that the return type of GetMomentumDirectionChange is a
  // pointer to G4ParticleMometum. Also it is a normalized 
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = theMassChange;

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
 
  // calculate new kinetic energy
  G4double energy = pPostStepPoint->GetKineticEnergy() 
                    + (theEnergyChange - pPreStepPoint->GetKineticEnergy()); 

  // update kinetic energy and momentum direction
  if (energy > 0.0) {
    // calculate new momentum
    G4ThreeVector pMomentum =  pPostStepPoint->GetMomentum() 
                + ( CalcMomentum(theEnergyChange, theMomentumDirectionChange, mass)
	            - pPreStepPoint->GetMomentum());
    G4double      tMomentum = pMomentum.mag();
    G4ThreeVector direction( pMomentum.x()/tMomentum,
			     pMomentum.y()/tMomentum,
			     pMomentum.z()/tMomentum );
    pPostStepPoint->SetMomentumDirection(direction);
    pPostStepPoint->SetKineticEnergy( energy );
  } else {
    // stop case
    pPostStepPoint->SetMomentumDirection(G4ThreeVector(1., 0., 0.));
    pPostStepPoint->SetKineticEnergy(0.0);
  }

  // update polarization
  pPostStepPoint->AddPolarization( thePolarizationChange
				   - pPreStepPoint->GetPolarization());
      
  // update position and time
  pPostStepPoint->AddPosition( thePositionChange 
			       - pPreStepPoint->GetPosition() );
  pPostStepPoint->AddGlobalTime( theTimeChange
				 - pPreStepPoint->GetGlobalTime());
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - pPreStepPoint->GetGlobalTime()); 
  pPostStepPoint->AddProperTime( theProperTimeChange 
				 - pPreStepPoint->GetProperTime());

  // update weight if use EB
  pPostStepPoint->SetWeight( theWeightChange );

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

G4Step* G4ParticleChange::UpdateStepForPostStep(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  // Take note that the return type of GetMomentumChange is a
  // pointer to G4ParticleMometum. Also it is a normalized 
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = theMassChange;

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

   // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - aTrack->GetGlobalTime());
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  // update weight if use EB
  pPostStepPoint->SetWeight( theWeightChange );

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}


G4Step* G4ParticleChange::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = theMassChange;

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - aTrack->GetGlobalTime());
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  // update weight if use EB
  pPostStepPoint->SetWeight( theWeightChange );

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

//----------------------------------------------------------------
// methods for printing messages  
//

void G4ParticleChange::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4cout.precision(3);

  G4cout << "        Mass (GeV)   : " 
       << setw(20) << theMassChange/GeV
       << endl; 
  G4cout << "        Charge (eplus)   : " 
       << setw(20) << theChargeChange/eplus
       << endl; 
  G4cout << "        Position - x (mm)   : " 
       << setw(20) << thePositionChange.x()/mm
       << endl; 
  G4cout << "        Position - y (mm)   : " 
       << setw(20) << thePositionChange.y()/mm
       << endl; 
  G4cout << "        Position - z (mm)   : " 
       << setw(20) << thePositionChange.z()/mm
       << endl;
  G4cout << "        Time (ns)           : " 
       << setw(20) << theTimeChange/ns
       << endl;
  G4cout << "        Proper Time (ns)    : " 
       << setw(20) << theProperTimeChange/ns
       << endl;
  G4cout << "        Momentum Direct - x : " 
       << setw(20) << theMomentumDirectionChange.x()
       << endl;
  G4cout << "        Momentum Direct - y : " 
       << setw(20) << theMomentumDirectionChange.y()
       << endl;
  G4cout << "        Momentum Direct - z : " 
       << setw(20) << theMomentumDirectionChange.z()
       << endl;
  G4cout << "        Kinetic Energy (MeV): " 
       << setw(20) << theEnergyChange/MeV
       << endl;
  G4cout << "        Polarization - x    : " 
       << setw(20) << thePolarizationChange.x()
       << endl;
  G4cout << "        Polarization - y    : " 
       << setw(20) << thePolarizationChange.y()
       << endl;
  G4cout << "        Polarization - z    : " 
       << setw(20) <<  thePolarizationChange.z()
       << endl;
  if (fUseEB) {
    G4cout << "        Track Weight      : " 
         << setw(20) <<  theWeightChange
         << endl;	
  }
}

G4bool G4ParticleChange::CheckIt(const G4Track& aTrack)
{
  G4bool    exitWithError = false;
  G4double  accuracy;

  // No check in case of "fStopAndKill" 
  if (GetStatusChange() ==   fStopAndKill )  {
    return G4VParticleChange::CheckIt(aTrack);
  }

  // MomentumDirection should be unit vector
  G4bool itsOKforMomentum = true;  
  if ( theEnergyChange >0.) {
    accuracy = abs(theMomentumDirectionChange.mag2()-1.0);
    if (accuracy > accuracyForWarning) {
      G4cout << "  G4ParticleChange::CheckIt  : ";
      G4cout << "the Momentum Change is not unit vector !!" << endl;
      G4cout << "  Difference:  " << accuracy << endl;
      itsOKforMomentum = false;
      if (accuracy > accuracyForException) exitWithError = true;
    }
  }

  // Both global and proper time should not go back
  G4bool itsOKforGlobalTime = true;  
  accuracy = (aTrack.GetGlobalTime()- theTimeChange)/ns;
  if (accuracy > accuracyForWarning) {
    G4cout << "  G4ParticleChange::CheckIt    : ";
    G4cout << "the global time goes back  !!" << endl;
    G4cout << "  Difference:  " << accuracy  << "[ns] " <<endl;
    itsOKforGlobalTime = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  G4bool itsOKforProperTime = true;
  accuracy = (aTrack.GetProperTime() - theProperTimeChange )/ns;
  if (accuracy > accuracyForWarning) {
    G4cout << "  G4ParticleChange::CheckIt    : ";
    G4cout << "the proper time goes back  !!" << endl;
    G4cout << "  Difference:  " << accuracy  << "[ns] " <<endl;
    itsOKforProperTime = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // Kinetic Energy should not be negative
  G4bool itsOKforEnergy = true;
  accuracy = -1.0*theEnergyChange/MeV;
  if (accuracy > accuracyForWarning) {
    G4cout << "  G4ParticleChange::CheckIt    : ";
    G4cout << "the kinetic energy is negative  !!" << endl;
    G4cout << "  Difference:  " << accuracy  << "[MeV] " <<endl;
    itsOKforEnergy = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  G4bool itsOK = itsOKforMomentum && itsOKforEnergy && itsOKforProperTime && itsOKforGlobalTime;
  // dump out information of this particle change
  if (!itsOK) { 
    G4cout << " G4ParticleChange::CheckIt " <<endl;
    G4cout << " pointer : " << this <<endl ;
    DumpInfo();
  }

  // Exit with error
  if (exitWithError) G4Exception("G4ParticleChange::CheckIt");

  //correction
  if (!itsOKforMomentum) {
    G4double vmag = theMomentumDirectionChange.mag();
    theMomentumDirectionChange = (1./vmag)*theMomentumDirectionChange;
  }
  if (!itsOKforGlobalTime) {
    theTimeChange = aTrack.GetGlobalTime();
  }
  if (!itsOKforProperTime) {
    theProperTimeChange = aTrack.GetProperTime();
  }
  if (!itsOKforEnergy) {
    theEnergyChange = 0.0;
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}



