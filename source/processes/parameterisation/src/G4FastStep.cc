// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastStep.cc,v 1.2 1999-04-14 14:25:35 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//$Id:
//---------------------------------------------------------------
//
//  G4FastStep.cc
//
//  Description:
//    Encapsulates a G4ParticleChange and insure friendly interface
//    methods to manage the primary/secondaries final state for 
//    Fast Simulation Models.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//    Apr 98: MoraDeFreitas - G4FastStep becomes the G4ParticleChange
//                      for the Fast Simulation Process.
//
//---------------------------------------------------------------

#include "G4FastStep.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

void G4FastStep::Initialize(const G4FastTrack& fastTrack)
{
  // keeps the fastTrack reference
  fFastTrack=&fastTrack;

  // currentTrack will be used to Initialize the other data members
  const G4Track& currentTrack = *(fFastTrack->GetPrimaryTrack());

  // use base class's method at first
  G4VParticleChange::Initialize(currentTrack);

  // set Energy/Momentum etc. equal to those of the parent particle
  const G4DynamicParticle*  pParticle = currentTrack.GetDynamicParticle();
  theEnergyChange          = pParticle->GetKineticEnergy();
  theMomentumChange        = pParticle->GetMomentumDirection();
  thePolarizationChange    = pParticle->GetPolarization();
  theProperTimeChange      = pParticle->GetProperTime();

  // set Position/Time etc. equal to those of the parent track
  thePositionChange      = currentTrack.GetPosition();
  theTimeChange          = currentTrack.GetGlobalTime();

  // switch off stepping hit invokation by default:
  theSteppingControlFlag = AvoidHitInvocation;

  // event biasing weigth:
  theWeightChange        = currentTrack.GetWeight();
}  

//----------------------------------------
// -- Set the StopAndKilled signal
// -- and put kinetic energy to 0.0. in the
// -- G4ParticleChange.
//----------------------------------------
void G4FastStep::KillPrimaryTrack()
{
  SetPrimaryTrackFinalKineticEnergy(0.) ;
  SetStatusChange(fStopAndKill) ;
}

//--------------------
//
//--------------------
void 
G4FastStep::
SetPrimaryTrackFinalPosition(const G4ThreeVector &position,
			     G4bool localCoordinates)
{
  // Compute the position coordinate in global
  // reference system if needed ...
  G4ThreeVector globalPosition = position;
  if (localCoordinates) 
    globalPosition = fFastTrack->GetInverseAffineTransformation()->
      TransformPoint(position);
  // ...and feed the globalPosition:
  thePositionChange = globalPosition;
}

//--------------------
//
//--------------------
void 
G4FastStep::
SetPrimaryTrackFinalMomentum(const G4ThreeVector &momentum,
			     G4bool localCoordinates)
{
  // Compute the momentum in global reference
  // system if needed ...
  G4ThreeVector globalMomentum = momentum;
  if (localCoordinates)
    globalMomentum = fFastTrack->GetInverseAffineTransformation()->
      TransformAxis(momentum);
  // ...and feed the globalMomentum (ensuring unitarity)
  SetMomentumChange(globalMomentum.unit());
}


//--------------------
//
//--------------------
void 
G4FastStep::
SetPrimaryTrackFinalKineticEnergyAndDirection(G4double kineticEnergy,
					      const G4ThreeVector &direction,
					      G4bool localCoordinates)
{
  // Compute global direction if needed...
  G4ThreeVector globalDirection = direction;
  if (localCoordinates)
    globalDirection =fFastTrack->GetInverseAffineTransformation()-> 
      TransformAxis(direction);
  // ...and feed the globalMomentum (ensuring unitarity)
  SetMomentumChange(globalDirection.unit());
  SetPrimaryTrackFinalKineticEnergy(kineticEnergy);
}

//--------------------
//
//--------------------
void 
G4FastStep::
SetPrimaryTrackFinalPolarization(const G4ThreeVector &polarization,
				 G4bool localCoordinates)
{
  // Compute polarization in global system if needed:
  G4ThreeVector globalPolarization(polarization);
  if (localCoordinates)
    globalPolarization = fFastTrack->GetInverseAffineTransformation()->
      TransformAxis(globalPolarization);  
  // Feed the particle globalPolarization:
  thePolarizationChange = globalPolarization;
}

//--------------------
//
//--------------------
G4Track* G4FastStep::
CreateSecondaryTrack(const G4DynamicParticle& dynamics,
		     G4ThreeVector polarization,
		     G4ThreeVector position,
		     G4double time,
		     G4bool localCoordinates     )
{
  G4DynamicParticle dummyDynamics(dynamics);
  
  // ------------------------------------------
  // Add the polarization to the dummyDynamics:
  // ------------------------------------------
  dummyDynamics.SetPolarization(polarization.x(),
				polarization.y(),
				polarization.z());
  
  return CreateSecondaryTrack(dummyDynamics, position, time, localCoordinates);
}

//--------------------
//
//--------------------
G4Track* G4FastStep::
CreateSecondaryTrack(const G4DynamicParticle& dynamics,
		     G4ThreeVector position,
		     G4double time,
		     G4bool localCoordinates     )
{
  // ----------------------------------------
  // Quantities in global coordinates system.
  //  
  // The allocated globalDynamics is deleted
  // by the destructor of the G4Track.
  // ----------------------------------------
  G4DynamicParticle* globalDynamics =
    new G4DynamicParticle(dynamics);
  G4ThreeVector globalPosition(position);
  
  // -----------------------------------
  // Convert to global system if needed:
  // -----------------------------------
  if (localCoordinates)
    {
      // -- Momentum Direction:
      globalDynamics->SetMomentumDirection(fFastTrack->
					   GetInverseAffineTransformation()->
					   TransformAxis(globalDynamics->
							 GetMomentumDirection()));
      // -- Polarization:
      G4ThreeVector globalPolarization;
      globalPolarization = fFastTrack->GetInverseAffineTransformation()->
	TransformAxis(globalDynamics->GetPolarization());
      globalDynamics->SetPolarization(
				      globalPolarization.x(),
				      globalPolarization.y(),
				      globalPolarization.z()
				      );
      
      // -- Position:
      globalPosition = fFastTrack->GetInverseAffineTransformation()->
	TransformPoint(globalPosition);
    }
  
  //-------------------------------------
  // Create the G4Track of the secondary:
  //-------------------------------------
  G4Track* secondary = new G4Track(
				   globalDynamics,
				   time,
				   globalPosition
				   );
  
  //-------------------------------
  // and feed the changes:
  //-------------------------------
  AddSecondary(secondary);
  
  //--------------------------------------
  // returns the pointer on the secondary:
  //--------------------------------------
  return secondary;
}

// G4FastStep should never be Initialized in this way
// but we must define it to avoid warnings.
void G4FastStep::Initialize(const G4Track&) {
  G4Exception("G4FastStep::Initialize(const G4Track&) should never be called,\nyou must use instead the G4FastStep::Initialize(const G4FastTrack&)\nmethod!");
}

G4FastStep::G4FastStep():G4VParticleChange()
{
  if (verboseLevel>2) {
    G4cerr << "G4FastStep::G4FastStep() " << endl;
  }
}

G4FastStep::~G4FastStep() 
{
  if (verboseLevel>2) {
    G4cerr << "G4FastStep::~G4FastStep() " << endl;
  }
}

// copy and assignment operators are implemented as "shallow copy"
G4FastStep::G4FastStep(const G4FastStep &right)
{
   *this = right;
}


G4FastStep & G4FastStep::operator=(const G4FastStep &right)
{
   if (this != &right)
   {
     G4VParticleChange::operator=(right);
     theListOfSecondaries          = right.theListOfSecondaries;
     theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
     theNumberOfSecondaries        = right.theNumberOfSecondaries;
     theStatusChange               = right.theStatusChange;
     theMomentumChange             = right.theMomentumChange;
     thePolarizationChange         = right.thePolarizationChange;
     thePositionChange             = right.thePositionChange;
     theTimeChange                 = right.theTimeChange;
     theEnergyChange               = right.theEnergyChange;
     theTrueStepLength             = right.theTrueStepLength;
     theLocalEnergyDeposit         = right.theLocalEnergyDeposit;
     theSteppingControlFlag        = right.theSteppingControlFlag;
     theWeightChange               = right.theWeightChange;
   }
   return *this;
}





G4bool G4FastStep::operator==(const G4FastStep &right) const
{
   return ((G4VParticleChange *)this == (G4VParticleChange *) &right);
}

G4bool G4FastStep::operator!=(const G4FastStep &right) const
{
   return ((G4VParticleChange *)this != (G4VParticleChange *) &right);
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4FastStep::UpdateStepForPostStep(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  // Take note that the return type of GetMomentumChange is a
  // pointer to G4ParticleMometum. Also it is a normalized 
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = aTrack->GetDynamicParticle()->GetMass();
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

   // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - aTrack->GetGlobalTime());
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  // update weight
  pPostStepPoint->SetWeight( theWeightChange );

  //#define FAST_STEP_DEBUG 1
#ifdef FAST_STEP_DEBUG
  CheckIt(*aTrack);
#endif
  
  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

G4Step* G4FastStep::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = mass = aTrack->GetDynamicParticle()->GetMass();
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - aTrack->GetGlobalTime());
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  // update weight
  pPostStepPoint->SetWeight( theWeightChange );

#ifdef FAST_STEP_DEBUG
  CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

//----------------------------------------------------------------
// methods for printing messages  
//

void G4FastStep::DumpInfo() const
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
  G4cout << "        Time (ns)           : " 
       << setw(20) << theTimeChange/ns
       << endl;
  G4cout << "        Proper Time (ns)    : " 
       << setw(20) << theProperTimeChange/ns
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
}

#ifdef FAST_STEP_DEBUG
void G4FastStep::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  if (theEnergyChange > aTrack.GetKineticEnergy()) {
    G4cout << " !!! the energy becomes larger than the initial energy !!!"
	   << " :  " << (theEnergyChange -aTrack.GetKineticEnergy())/MeV
	   << "MeV " <<endl;
    itsOK = false;
  }
  if ( (theEnergyChange >0.) && 
       ( abs(theMomentumChange.mag2()-1.0) > 1.0e-5 ) ){
    G4cout << " !!! the Momentum Change is not unit vector !!!!"
	   << " :  " << theMomentumChange.mag()
	   << endl;
    itsOK = false;
  }
  if (theTimeChange < aTrack.GetGlobalTime()) {
    G4cout << " !!! the global time goes back  !!!"
	   << " :  " << aTrack.GetGlobalTime()/ns
	   << " -> " << theTimeChange/ns
	   << "[ns] " <<endl;
    itsOK = false;
  }
  if (theProperTimeChange < aTrack.GetProperTime()) {
    G4cout << " !!! the poper time goes back  !!!"
	   << " :  " << aTrack.GetProperTime()/ns
	   << " -> " << theProperTimeChange/ns
	   << "[ns] " <<endl;
    itsOK = false;
  }
  
  if (!itsOK) { 
    G4cout << " G4FastStep::CheckIt " <<endl;
    G4cout << " pointer : " << this <<endl ;
    DumpInfo();
    G4Exception("\nG4FastStep was compiled with FAST_STEP_DEBUG option.\nG4FastStep::CheckIt detected a fatal error. \nPlease report the problem via the Geant4 Web site,\naddress http://wwwinfo.cern.ch/asd/geant4.");
  }
}

#endif
