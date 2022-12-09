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
//
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

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
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
  ProposeTrackStatus(fStopAndKill) ;
}

//--------------------
//
//--------------------
void 
G4FastStep::
ProposePrimaryTrackFinalPosition(const G4ThreeVector &position,
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

void 
G4FastStep::
SetPrimaryTrackFinalPosition(const G4ThreeVector &position,
                             G4bool localCoordinates)
{
  ProposePrimaryTrackFinalPosition(position, localCoordinates);
}

//--------------------
//
//--------------------
void 
G4FastStep::
ProposePrimaryTrackFinalMomentumDirection(const G4ThreeVector &momentum,
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

void 
G4FastStep::
SetPrimaryTrackFinalMomentum(const G4ThreeVector &momentum,
                             G4bool localCoordinates)
{
  ProposePrimaryTrackFinalMomentumDirection(momentum, localCoordinates);
}

//--------------------
//
//--------------------
void 
G4FastStep::
ProposePrimaryTrackFinalKineticEnergyAndDirection(G4double kineticEnergy,
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

void 
G4FastStep::
SetPrimaryTrackFinalKineticEnergyAndDirection(G4double kineticEnergy,
                                              const G4ThreeVector &direction,
                                              G4bool localCoordinates)
{
  ProposePrimaryTrackFinalKineticEnergyAndDirection(kineticEnergy, direction, localCoordinates);
}

//--------------------
//
//--------------------
void 
G4FastStep::
ProposePrimaryTrackFinalPolarization(const G4ThreeVector &polarization,
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

void 
G4FastStep::
SetPrimaryTrackFinalPolarization(const G4ThreeVector &polarization,
                                 G4bool localCoordinates)
{
  ProposePrimaryTrackFinalPolarization(polarization, localCoordinates);
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
void G4FastStep::Initialize(const G4Track&)
{
  G4ExceptionDescription  tellWhatIsWrong;
  tellWhatIsWrong << "G4FastStep can be initialised only through G4FastTrack."
		  << G4endl;
  G4Exception("G4FastStep::Initialize(const G4Track&)",
	      "FastSim005",
	      FatalException,
	      tellWhatIsWrong);
}

G4FastStep::G4FastStep()
  : G4VParticleChange()
{
  if (verboseLevel>2)
  {
    G4cerr << "G4FastStep::G4FastStep()" << G4endl;
  }
}

G4FastStep::~G4FastStep() 
{
  if (verboseLevel>2)
  {
    G4cerr << "G4FastStep::~G4FastStep()" << G4endl;
  }
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

  //  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  //  G4double     mass = aTrack->GetDynamicParticle()->GetMass();
 
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

  if (debugFlag) CheckIt(*aTrack);

  
  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

G4Step* G4FastStep::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  // G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();
  // G4double     mass = aTrack->GetDynamicParticle()->GetMass();
 
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

  if (debugFlag) CheckIt(*aTrack);

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
  
  G4cout << "        Position - x (mm)   : " <<      G4BestUnit( thePositionChange.x(), "Length" ) << G4endl; 
  G4cout << "        Position - y (mm)   : " <<      G4BestUnit( thePositionChange.y(), "Length" ) << G4endl;
  G4cout << "        Position - z (mm)   : " <<      G4BestUnit( thePositionChange.z(), "Length" ) << G4endl;   
  G4cout << "        Time (ns)           : " <<      G4BestUnit( theTimeChange,         "Time"   ) << G4endl;
  G4cout << "        Proper Time (ns)    : " <<      G4BestUnit( theProperTimeChange,   "Time"   ) << G4endl;
  G4long olprc = G4cout.precision(3);
  G4cout << "        Momentum Direct - x : " << std::setw(20) << theMomentumChange.x() << G4endl;
  G4cout << "        Momentum Direct - y : " << std::setw(20) << theMomentumChange.y() << G4endl;
  G4cout << "        Momentum Direct - z : " << std::setw(20) << theMomentumChange.z() << G4endl;
  G4cout.precision(olprc);
  G4cout << "        Kinetic Energy (MeV): " <<      G4BestUnit( theEnergyChange,       "Energy" ) << G4endl;
  G4cout.precision(3);
  G4cout << "        Polarization - x    : " << std::setw(20) << thePolarizationChange.x() << G4endl;
  G4cout << "        Polarization - y    : " << std::setw(20) << thePolarizationChange.y() << G4endl;
  G4cout << "        Polarization - z    : " << std::setw(20) << thePolarizationChange.z() << G4endl;
  G4cout.precision(olprc);
}

G4bool G4FastStep::CheckIt(const G4Track& aTrack)
{
  //
  //      In the G4FastStep::CheckIt
  //      We only check a bit
  //      
  //      If the user violates the energy,
  //      We don't care, we agree.
  //
  //      But for theMomentumDirectionChange,
  //      We do pay attention.
  //      And if too large is its range,
  //      We issue an Exception.
  //
  //
  // It means, the G4FastStep::CheckIt issues an exception only for the
  // theMomentumDirectionChange which should be an unit vector
  // and it corrects it because it could cause problems for the ulterior
  // tracking.For the rest, only warning are issued.

  G4bool    itsOK = true;
  G4bool    exitWithError = false;
  G4double  accuracy;
  
  // Energy should not be larger than the initial value
  accuracy = ( theEnergyChange - aTrack.GetKineticEnergy())/MeV;
  if (accuracy > GetAccuracyForWarning())
    {
      G4ExceptionDescription ed;
      ed << "The energy becomes larger than the initial value, difference = " <<  accuracy << " MeV" << G4endl;
      G4Exception("G4FastStep::CheckIt(const G4Track& aTrack)",
		  "FastSim006",
		  JustWarning, ed);
      itsOK = false;
      if (accuracy > GetAccuracyForException())  {exitWithError = true;}
    }
  
  G4bool itsOKforMomentum = true;
  if ( theEnergyChange >0.)
    {
      accuracy = std::abs(theMomentumChange.mag2()-1.0);
      if (accuracy > GetAccuracyForWarning())
	{
	  G4ExceptionDescription ed;
	  ed << "The Momentum Change is not a unit vector, difference = " <<  accuracy << G4endl;
	  G4Exception("G4FastStep::CheckIt(const G4Track& aTrack)",
		      "FastSim007",
		      JustWarning, ed);
	  itsOK = itsOKforMomentum = false;
	  if (accuracy > GetAccuracyForException())  {exitWithError = true;}
	}
    }
  
  accuracy = (aTrack.GetGlobalTime()- theTimeChange)/ns;  
  if (accuracy > GetAccuracyForWarning())
    {
      G4ExceptionDescription ed;
      ed << "The global time is getting backward, difference = " <<  accuracy << " ns" << G4endl;
      G4Exception("G4FastStep::CheckIt(const G4Track& aTrack)",
		  "FastSim008",
		  JustWarning, ed);
      itsOK = false;
    }
  
  accuracy = (aTrack.GetProperTime() - theProperTimeChange )/ns;
  if (accuracy >  GetAccuracyForWarning())
    {
      G4ExceptionDescription ed;
      ed << "The proper time is getting backward, difference = " <<  accuracy << " ns" << G4endl;
      G4Exception("G4FastStep::CheckIt(const G4Track& aTrack)",
		  "FastSim009",
		  JustWarning, ed);
      itsOK = false;
    }
  
  if (!itsOK)
    { 
      G4cout << "ERROR - G4FastStep::CheckIt() " << G4endl;
      G4cout << "        Pointer : " << this << G4endl ;
      DumpInfo();
    }
  
  // Exit with error
  if (exitWithError)
    {
      G4ExceptionDescription ed;
      ed << "An inaccuracy in G4FastStep is beyond tolerance." << G4endl;
      G4Exception("G4FastStep::CheckIt(const G4Track& aTrack)",
		  "FastSim010",
		  FatalException, ed);
    }
  
  //correction for Momentum only.
  if (!itsOKforMomentum) {
    G4double vmag = theMomentumChange.mag();
    theMomentumChange = (1./vmag)*theMomentumChange;
  }
  
  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack); 
  return itsOK;
}
