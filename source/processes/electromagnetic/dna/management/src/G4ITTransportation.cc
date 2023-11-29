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
/// \brief This class is a slightly modified version of G4Transportation
///        initially written by John Apostolakis and colleagues
///        But it should use the exact same algorithm
//
// Contact : Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History :
// -----------
// =======================================================================
// Modified:
//   28 Oct  2011, P.Gumpl./J.Ap: Detect gravity field, use magnetic moment
//   20 Nov  2008, J.Apostolakis: Push safety to helper - after ComputeSafety
//    9 Nov  2007, J.Apostolakis: Flag for short steps, push safety to helper
//   19 Jan  2006, P.MoraDeFreitas: Fix for suspended tracks (StartTracking)
//   11 Aug  2004, M.Asai: Add G4VSensitiveDetector* for updating stepPoint.
//   21 June 2003, J.Apostolakis: Calling field manager with
//                     track, to enable it to configure its accuracy
//   13 May  2003, J.Apostolakis: Zero field areas now taken into
//                     account correclty in all cases (thanks to W Pokorski).
//   29 June 2001, J.Apostolakis, D.Cote-Ahern, P.Gumplinger:
//                     correction for spin tracking
//   20 Febr 2001, J.Apostolakis:  update for new FieldTrack
//   22 Sept 2000, V.Grichine:     update of Kinetic Energy
//  ---------------------------------------------------
//   10 Oct  2011, M.Karamitros:   G4ITTransportation created
// Created:  19 March 1997, J. Apostolakis
// =======================================================================
//
// -------------------------------------------------------------------

#include "G4ITTransportation.hh"
#include "G4IT.hh"
#include "G4TrackingInformation.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4ITTransportationManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ITNavigator.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4ITSafetyHelper.hh"
#include "G4FieldManagerStore.hh"
#include "G4LowEnergyEmProcessSubType.hh"

#include "G4UnitsTable.hh"
#include "G4ReferenceCast.hh"

class G4VSensitiveDetector;

#ifndef PrepareState
#  define PrepareState()                                                       \
    G4ITTransportationState* __state = this->GetState<G4ITTransportationState>()
#endif

#ifndef State
#define State(theXInfo) (__state->theXInfo)
#endif

//#define DEBUG_MEM

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
using G4MemStat::MemStat;
#endif

//#define G4DEBUG_TRANSPORT 1

G4ITTransportation::G4ITTransportation(const G4String& aName, int verbose) :
    G4VITProcess(aName, fTransportation),
    fThreshold_Warning_Energy(100 * MeV),
    fThreshold_Important_Energy(250 * MeV),
    fThresholdTrials(10),
    fUnimportant_Energy(1 * MeV), // Not used
    fSumEnergyKilled(0.0),
    fMaxEnergyKilled(0.0),
    fShortStepOptimisation(false), // Old default: true (=fast short steps)
    fVerboseLevel(verbose)
{
  pParticleChange = &fParticleChange;
  G4TransportationManager* transportMgr;
  transportMgr = G4TransportationManager::GetTransportationManager();
  G4ITTransportationManager* ITtransportMgr;
  ITtransportMgr = G4ITTransportationManager::GetTransportationManager();
  fLinearNavigator = ITtransportMgr->GetNavigatorForTracking();
  fFieldPropagator = transportMgr->GetPropagatorInField();
  fpSafetyHelper = ITtransportMgr->GetSafetyHelper(); // New

  // Cannot determine whether a field exists here, as it would
  //  depend on the relative order of creating the detector's
  //  field and this process. That order is not guaranted.
  // Instead later the method DoesGlobalFieldExist() is called

  enableAtRestDoIt = false;
  enableAlongStepDoIt = true;
  enablePostStepDoIt = true;
  SetProcessSubType(fLowEnergyTransportation);
  SetInstantiateProcessState(true);
  G4VITProcess::SetInstantiateProcessState(false);
  fInstantiateProcessState = true;

  G4VITProcess::fpState.reset(new G4ITTransportationState());
  /*
   if(fTransportationState == 0)
   {
   G4cout << "KILL in G4ITTransportation::G4ITTransportation" << G4endl;
   abort();
   }
   */
}

G4ITTransportation::G4ITTransportation(const G4ITTransportation& right) :
    G4VITProcess(right)
{
  // Copy attributes
  fVerboseLevel = right.fVerboseLevel;
  fThreshold_Warning_Energy = right.fThreshold_Warning_Energy;
  fThreshold_Important_Energy = right.fThreshold_Important_Energy;
  fThresholdTrials = right.fThresholdTrials;
  fUnimportant_Energy = right.fUnimportant_Energy;
  fSumEnergyKilled = right.fSumEnergyKilled;
  fMaxEnergyKilled = right.fMaxEnergyKilled;
  fShortStepOptimisation = right.fShortStepOptimisation;

  // Setup Navigators
  G4TransportationManager* transportMgr;
  transportMgr = G4TransportationManager::GetTransportationManager();
  G4ITTransportationManager* ITtransportMgr;
  ITtransportMgr = G4ITTransportationManager::GetTransportationManager();
  fLinearNavigator = ITtransportMgr->GetNavigatorForTracking();
  fFieldPropagator = transportMgr->GetPropagatorInField();
  fpSafetyHelper = ITtransportMgr->GetSafetyHelper(); // New

  // Cannot determine whether a field exists here, as it would
  //  depend on the relative order of creating the detector's
  //  field and this process. That order is not guaranted.
  // Instead later the method DoesGlobalFieldExist() is called

  enableAtRestDoIt = false;
  enableAlongStepDoIt = true;
  enablePostStepDoIt = true;

  pParticleChange = &fParticleChange;
  SetInstantiateProcessState(true);
  G4VITProcess::SetInstantiateProcessState(false);
  fInstantiateProcessState = right.fInstantiateProcessState;
}

G4ITTransportation& G4ITTransportation::operator=(const G4ITTransportation& /*right*/)
{
//  if (this == &right) return *this;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////
/// Process State
//////////////////////////////////////////////////////////////////////////////
G4ITTransportation::G4ITTransportationState::G4ITTransportationState() :
    G4ProcessState(), fCurrentTouchableHandle(0)
{
  fTransportEndPosition = G4ThreeVector(0, 0, 0);
  fTransportEndMomentumDir = G4ThreeVector(0, 0, 0);
  fTransportEndKineticEnergy = -1;
  fTransportEndSpin = G4ThreeVector(0, 0, 0);
  fMomentumChanged = false;
  fEnergyChanged = false;
  fEndGlobalTimeComputed = false;
  fCandidateEndGlobalTime = -1;
  fParticleIsLooping = false;

  static G4ThreadLocal G4TouchableHandle *nullTouchableHandle = 0;
  if (!nullTouchableHandle) nullTouchableHandle = new G4TouchableHandle;
  // Points to (G4VTouchable*) 0

  fCurrentTouchableHandle = *nullTouchableHandle;
  fGeometryLimitedStep = false;
  fPreviousSftOrigin = G4ThreeVector(0, 0, 0);
  fPreviousSafety = 0.0;
  fNoLooperTrials = false;
  fEndPointDistance = -1;
}

G4ITTransportation::G4ITTransportationState::~G4ITTransportationState()
{
  ;
}

G4ITTransportation::~G4ITTransportation()
{
#ifdef G4VERBOSE
  if ((fVerboseLevel > 0) && (fSumEnergyKilled > 0.0))
  {
    G4cout << " G4ITTransportation: Statistics for looping particles "
           << G4endl;
    G4cout << "   Sum of energy of loopers killed: "
           << fSumEnergyKilled << G4endl;
    G4cout << "   Max energy of loopers killed: "
           << fMaxEnergyKilled << G4endl;
  }
#endif
}

void G4ITTransportation::BuildPhysicsTable(const G4ParticleDefinition&)
{
  fpSafetyHelper->InitialiseHelper();
}

G4bool G4ITTransportation::DoesGlobalFieldExist()
{
  G4TransportationManager* transportMgr;
  transportMgr = G4TransportationManager::GetTransportationManager();

  // fFieldExists= transportMgr->GetFieldManager()->DoesFieldExist();
  // return fFieldExists;
  return transportMgr->GetFieldManager()->DoesFieldExist();
}

//////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.
G4double
G4ITTransportation::
AlongStepGetPhysicalInteractionLength(const G4Track& track,
                                      G4double,
                                      G4double currentMinimumStep,
                                      G4double& currentSafety,
                                      G4GPILSelection* selection)
{
  PrepareState();
  G4double geometryStepLength(-1.0), newSafety(-1.0);

  State(fParticleIsLooping) = false;
  State(fEndGlobalTimeComputed) = false;
  State(fGeometryLimitedStep) = false;

  // Initial actions moved to  StartTrack()
  // --------------------------------------
  // Note: in case another process changes touchable handle
  //    it will be necessary to add here (for all steps)
  // State(fCurrentTouchableHandle) = track.GetTouchableHandle();

  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  //
  *selection = CandidateForSelection;

  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle* pParticle = track.GetDynamicParticle();
//    const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
  G4ThreeVector startMomentumDir = pParticle->GetMomentumDirection();
  G4ThreeVector startPosition = track.GetPosition();

  // G4double   theTime        = track.GetGlobalTime() ;

  // The Step Point safety can be limited by other geometries and/or the
  // assumptions of any process - it's not always the geometrical safety.
  // We calculate the starting point's isotropic safety here.
  //
  G4ThreeVector OriginShift = startPosition - State(fPreviousSftOrigin);
  G4double MagSqShift = OriginShift.mag2();
  if (MagSqShift >= sqr(State(fPreviousSafety)))
  {
    currentSafety = 0.0;
  }
  else
  {
    currentSafety = State(fPreviousSafety) - std::sqrt(MagSqShift);
  }

  // Is the particle charged ?
  //
  G4double particleCharge = pParticle->GetCharge();

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction
  //   By the tracking, after all AlongStepDoIts, in "Relocation"

  // Check whether the particle have an (EM) field force exerting upon it
  //
  G4FieldManager* fieldMgr = 0;
  G4bool fieldExertsForce = false;
  if ((particleCharge != 0.0))
  {
    fieldMgr = fFieldPropagator->FindAndSetFieldManager(track.GetVolume());
    if (fieldMgr != 0)
    {
      // Message the field Manager, to configure it for this track
      fieldMgr->ConfigureForTrack(&track);
      // Moved here, in order to allow a transition
      //   from a zero-field  status (with fieldMgr->(field)0
      //   to a finite field  status

      // If the field manager has no field, there is no field !
      fieldExertsForce = (fieldMgr->GetDetectorField() != 0);
    }
  }

  // G4cout << " G4Transport:  field exerts force= " << fieldExertsForce
  // 	 << "  fieldMgr= " << fieldMgr << G4endl;

  // Choose the calculation of the transportation: Field or not
  //
  if (!fieldExertsForce)
  {
    G4double linearStepLength;
    if (fShortStepOptimisation && (currentMinimumStep <= currentSafety))
    {
      // The Step is guaranteed to be taken
      //
      geometryStepLength = currentMinimumStep;
      State(fGeometryLimitedStep) = false;
    }
    else
    {
      //  Find whether the straight path intersects a volume
      //
      // fLinearNavigator->SetNavigatorState(GetIT(track)->GetTrackingInfo()->GetNavigatorState());
      linearStepLength = fLinearNavigator->ComputeStep(startPosition,
                                                       startMomentumDir,
                                                       currentMinimumStep,
                                                       newSafety);

//      G4cout << "linearStepLength : " <<  G4BestUnit(linearStepLength,"Length")
//          << " | currentMinimumStep: " << currentMinimumStep
//          << " | trackID: " << track.GetTrackID() << G4endl;

      // Remember last safety origin & value.
      //
      State(fPreviousSftOrigin) = startPosition;
      State(fPreviousSafety) = newSafety;
      
      G4TrackStateManager& trackStateMan = GetIT(track)->GetTrackingInfo()
          ->GetTrackStateManager();
      fpSafetyHelper->LoadTrackState(trackStateMan);
      // fpSafetyHelper->SetTrackState(state);
      fpSafetyHelper->SetCurrentSafety(newSafety,
                                       State(fTransportEndPosition));
      fpSafetyHelper->ResetTrackState();
      
      // The safety at the initial point has been re-calculated:
      //
      currentSafety = newSafety;

      State(fGeometryLimitedStep) = (linearStepLength <= currentMinimumStep);
      if (State(fGeometryLimitedStep))
      {
        // The geometry limits the Step size (an intersection was found.)
        geometryStepLength = linearStepLength;
      }
      else
      {
        // The full Step is taken.
        geometryStepLength = currentMinimumStep;
      }
    }
    State(fEndPointDistance) = geometryStepLength;

    // Calculate final position
    //
    State(fTransportEndPosition) = startPosition
        + geometryStepLength * startMomentumDir;

    // Momentum direction, energy and polarisation are unchanged by transport
    //
    State(fTransportEndMomentumDir) = startMomentumDir;
    State(fTransportEndKineticEnergy) = track.GetKineticEnergy();
    State(fTransportEndSpin) = track.GetPolarization();
    State(fParticleIsLooping) = false;
    State(fMomentumChanged) = false;
    State(fEndGlobalTimeComputed) = true;
    State(theInteractionTimeLeft) = State(fEndPointDistance)
        / track.GetVelocity();
    State(fCandidateEndGlobalTime) = State(theInteractionTimeLeft)
        + track.GetGlobalTime();
/*
    G4cout << "track.GetVelocity() : "
           << track.GetVelocity() << G4endl;
    G4cout << "State(endpointDistance) : "
           << G4BestUnit(State(endpointDistance),"Length") << G4endl;
    G4cout << "State(theInteractionTimeLeft) : "
           << G4BestUnit(State(theInteractionTimeLeft),"Time") << G4endl;
    G4cout << "track.GetGlobalTime() : "
           << G4BestUnit(track.GetGlobalTime(),"Time") << G4endl;
*/
  }
  else //  A field exerts force
  {

    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "ITTransportation does not support external fields.";
    exceptionDescription
        << " If you are dealing with a tradiational MC simulation, ";
    exceptionDescription << "please use G4Transportation.";

    G4Exception("G4ITTransportation::AlongStepGetPhysicalInteractionLength",
                "NoExternalFieldSupport", FatalException, exceptionDescription);
    /*
     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
     //        G4ThreeVector  EndUnitMomentum ;
     G4double       lengthAlongCurve ;
     G4double       restMass = pParticleDef->GetPDGMass() ;

     fFieldPropagator->SetChargeMomentumMass( particleCharge,    // in e+ units
     momentumMagnitude, // in Mev/c
     restMass           ) ;

     G4ThreeVector spin        = track.GetPolarization() ;
     G4FieldTrack  aFieldTrack = G4FieldTrack( startPosition,
     track.GetMomentumDirection(),
     0.0,
     track.GetKineticEnergy(),
     restMass,
     track.GetVelocity(),
     track.GetGlobalTime(), // Lab.
     track.GetProperTime(), // Part.
     &spin                  ) ;
     if( currentMinimumStep > 0 )
     {
     // Do the Transport in the field (non recti-linear)
     //
     lengthAlongCurve = fFieldPropagator->ComputeStep( aFieldTrack,
     currentMinimumStep,
     currentSafety,
     track.GetVolume() ) ;
     State(fGeometryLimitedStep)= lengthAlongCurve < currentMinimumStep;
     if( State(fGeometryLimitedStep) )
     {
     geometryStepLength   = lengthAlongCurve ;
     }
     else
     {
     geometryStepLength   = currentMinimumStep ;
     }

     // Remember last safety origin & value.
     //
     State(fPreviousSftOrigin) = startPosition ;
     State(fPreviousSafety)    = currentSafety ;
     fpSafetyHelper->SetCurrentSafety( newSafety, startPosition);
     }
     else
     {
     geometryStepLength   = lengthAlongCurve= 0.0 ;
     State(fGeometryLimitedStep) = false ;
     }

     // Get the End-Position and End-Momentum (Dir-ection)
     //
     State(fTransportEndPosition) = aFieldTrack.GetPosition() ;

     // Momentum:  Magnitude and direction can be changed too now ...
     //
     State(fMomentumChanged)         = true ;
     State(fTransportEndMomentumDir) = aFieldTrack.GetMomentumDir() ;

     State(fTransportEndKineticEnergy)  = aFieldTrack.GetKineticEnergy() ;

     if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() )
     {
     // If the field can change energy, then the time must be integrated
     //    - so this should have been updated
     //
     State(fCandidateEndGlobalTime)   = aFieldTrack.GetLabTimeOfFlight();
     State(fEndGlobalTimeComputed)    = true;

     State(theInteractionTimeLeft) = State(fCandidateEndGlobalTime) -
                                         track.GetGlobalTime() ;

     // was ( State(fCandidateEndGlobalTime) != track.GetGlobalTime() );
     // a cleaner way is to have FieldTrack knowing whether time is updated.
     }
     else
     {
     // The energy should be unchanged by field transport,
     //    - so the time changed will be calculated elsewhere
     //
     State(fEndGlobalTimeComputed) = false;

     // Check that the integration preserved the energy
     //     -  and if not correct this!
     G4double  startEnergy= track.GetKineticEnergy();
     G4double  endEnergy= State(fTransportEndKineticEnergy);

     static G4int no_inexact_steps=0, no_large_ediff;
     G4double absEdiff = std::fabs(startEnergy- endEnergy);
     if( absEdiff > perMillion * endEnergy )
     {
     no_inexact_steps++;
     // Possible statistics keeping here ...
     }
     #ifdef G4VERBOSE
     if( fVerboseLevel > 1 )
     {
     if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
     {
     static G4int no_warnings= 0, warnModulo=1,  moduloFactor= 10;
     no_large_ediff ++;
     if( (no_large_ediff% warnModulo) == 0 )
     {
     no_warnings++;
     G4cout << "WARNING - G4Transportation::AlongStepGetPIL() "
     << "   Energy change in Step is above 1^-3 relative value. " << G4endl
     << "   Relative change in 'tracking' step = "
     << std::setw(15) << (endEnergy-startEnergy)/startEnergy << G4endl
     << "     Starting E= " << std::setw(12) << startEnergy / MeV << " MeV "
     << G4endl
     << "     Ending   E= " << std::setw(12) << endEnergy   / MeV << " MeV "
     << G4endl;
     G4cout << " Energy has been corrected -- however, review"
     << " field propagation parameters for accuracy."  << G4endl;
     if( (fVerboseLevel > 2 ) || (no_warnings<4) ||
         (no_large_ediff == warnModulo * moduloFactor) )
     {
     G4cout << " These include EpsilonStepMax(/Min) in G4FieldManager "
     << " which determine fractional error per step for integrated quantities. "
     << G4endl
     << " Note also the influence of the permitted number of integration steps."
     << G4endl;
     }
     G4cerr << "ERROR - G4Transportation::AlongStepGetPIL()" << G4endl
     << "        Bad 'endpoint'. Energy change detected"
     << " and corrected. "
     << " Has occurred already "
     << no_large_ediff << " times." << G4endl;
     if( no_large_ediff == warnModulo * moduloFactor )
     {
     warnModulo *= moduloFactor;
     }
     }
     }
     }  // end of if (fVerboseLevel)
     #endif
     // Correct the energy for fields that conserve it
     //  This - hides the integration error
     //       - but gives a better physical answer
     State(fTransportEndKineticEnergy)= track.GetKineticEnergy();
     }

     State(fTransportEndSpin) = aFieldTrack.GetSpin();
     State(fParticleIsLooping) = fFieldPropagator->IsParticleLooping() ;
     State(endpointDistance)   = (State(fTransportEndPosition) -
                                  startPosition).mag() ;
     // State(theInteractionTimeLeft) =
                       track.GetVelocity()/State(endpointDistance) ;
     */
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if (currentMinimumStep == 0.0)
  {
    if (currentSafety == 0.0)
    {
      State(fGeometryLimitedStep) = true;
//            G4cout << "!!!! Safety is NULL, on the Boundary !!!!!" << G4endl;
//            G4cout << " Track position : " << track.GetPosition() /nanometer
//                   << G4endl;
    }
  }

  // Update the safety starting from the end-point,
  // if it will become negative at the end-point.
  //
  if (currentSafety < State(fEndPointDistance))
  {
    // if( particleCharge == 0.0 )
    //    G4cout  << "  Avoiding call to ComputeSafety : charge = 0.0 " << G4endl;

    if (particleCharge != 0.0)
    {

      G4double endSafety = fLinearNavigator->ComputeSafety(
          State(fTransportEndPosition));
      currentSafety = endSafety;
      State(fPreviousSftOrigin) = State(fTransportEndPosition);
      State(fPreviousSafety) = currentSafety;
      
      /*
       G4VTrackStateHandle state =
       GetIT(track)->GetTrackingInfo()->GetTrackState(fpSafetyHelper);
       */
      G4TrackStateManager& trackStateMan = GetIT(track)->GetTrackingInfo()
          ->GetTrackStateManager();
      fpSafetyHelper->LoadTrackState(trackStateMan);
      // fpSafetyHelper->SetTrackState(state);
      fpSafetyHelper->SetCurrentSafety(currentSafety,
                                       State(fTransportEndPosition));
      fpSafetyHelper->ResetTrackState();

      // Because the Stepping Manager assumes it is from the start point,
      //  add the StepLength
      //
      currentSafety += State(fEndPointDistance);

#ifdef G4DEBUG_TRANSPORT
      G4cout.precision(12);
      G4cout << "***G4Transportation::AlongStepGPIL ** " << G4endl;
      G4cout << "  Called Navigator->ComputeSafety at "
          << State(fTransportEndPosition)
          << "    and it returned safety= " << endSafety << G4endl;
      G4cout << "  Adding endpoint distance " << State(fEndPointDistance)
             << "   to obtain pseudo-safety= " << currentSafety << G4endl;
#endif
    }
  }

  // fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

//  G4cout << "G4ITTransportation::AlongStepGetPhysicalInteractionLength = "
//         << G4BestUnit(geometryStepLength,"Length") << G4endl;

  return geometryStepLength;
}

void G4ITTransportation::ComputeStep(const G4Track& track,
                                     const G4Step& /*step*/,
                                     const double timeStep,
                                     double& oPhysicalStep)
{
  PrepareState();
  const G4DynamicParticle* pParticle = track.GetDynamicParticle();
  G4ThreeVector startMomentumDir = pParticle->GetMomentumDirection();
  G4ThreeVector startPosition = track.GetPosition();

  track.CalculateVelocity();
  G4double initialVelocity = track.GetVelocity();

  State(fGeometryLimitedStep) = false;

  /////////////////////////
  // !!! CASE NO FIELD !!!
  /////////////////////////
  State(fCandidateEndGlobalTime) = timeStep + track.GetGlobalTime();
  State(fEndGlobalTimeComputed) = true;

  // Choose the calculation of the transportation: Field or not
  //
  if (!State(fMomentumChanged))
  {
    //        G4cout << "Momentum has not changed" << G4endl;
    fParticleChange.ProposeVelocity(initialVelocity);
    oPhysicalStep = initialVelocity * timeStep;

    // Calculate final position
    //
    State(fTransportEndPosition) = startPosition
        + oPhysicalStep * startMomentumDir;
  }
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)
#include "G4ParticleTable.hh"
G4VParticleChange* G4ITTransportation::AlongStepDoIt(const G4Track& track,
                                                     const G4Step& stepData)
{

#if defined (DEBUG_MEM)
  MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM)
  mem_first = MemoryUsage();
#endif

  PrepareState();

  // G4cout << "G4ITTransportation::AlongStepDoIt" << G4endl;
  // set  pdefOpticalPhoton
  // Andrea Dotti: the following statement should be in a single line:
  // G4-MT transformation tools get confused if statement spans two lines
  // If needed contact: adotti@slac.stanford.edu
  static G4ThreadLocal G4ParticleDefinition* pdefOpticalPhoton = 0;
  if (!pdefOpticalPhoton) pdefOpticalPhoton =
      G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

  static G4ThreadLocal G4int noCalls = 0;
  noCalls++;

  fParticleChange.Initialize(track);

  //  Code for specific process
  //
  fParticleChange.ProposePosition(State(fTransportEndPosition));
  fParticleChange.ProposeMomentumDirection(State(fTransportEndMomentumDir));
  fParticleChange.ProposeEnergy(State(fTransportEndKineticEnergy));
  fParticleChange.SetMomentumChanged(State(fMomentumChanged));

  fParticleChange.ProposePolarization(State(fTransportEndSpin));

  G4double deltaTime = 0.0;

  // Calculate  Lab Time of Flight (ONLY if field Equations used it!)
  // G4double endTime   = State(fCandidateEndGlobalTime);
  // G4double delta_time = endTime - startTime;

  G4double startTime = track.GetGlobalTime();
  ///___________________________________________________________________________
  /// !!!!!!!
  /// A REVOIR !!!!
  if (State(fEndGlobalTimeComputed) == false)
  {
    // The time was not integrated .. make the best estimate possible
    //
    G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity();
    G4double stepLength = track.GetStepLength();

    deltaTime = 0.0; // in case initialVelocity = 0
    if (track.GetParticleDefinition() == pdefOpticalPhoton)
    {
      // For only Optical Photon, final velocity is used
      double finalVelocity = track.CalculateVelocityForOpticalPhoton();
      fParticleChange.ProposeVelocity(finalVelocity);
      deltaTime = stepLength / finalVelocity;
    }
    else if (initialVelocity > 0.0)
    {
      deltaTime = stepLength / initialVelocity;
    }

    State(fCandidateEndGlobalTime) = startTime + deltaTime;
  }
  else
  {
    deltaTime = State(fCandidateEndGlobalTime) - startTime;
  }

  fParticleChange.ProposeGlobalTime(State(fCandidateEndGlobalTime));
  fParticleChange.ProposeLocalTime(track.GetLocalTime() + deltaTime);
  /*
   // Now Correct by Lorentz factor to get delta "proper" Time

   G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
   G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;

   fParticleChange.ProposeProperTime(track.GetProperTime() + deltaProperTime) ;
   */

  fParticleChange.ProposeTrueStepLength(track.GetStepLength());

  ///___________________________________________________________________________
  ///

  // If the particle is caught looping or is stuck (in very difficult
  // boundaries) in a magnetic field (doing many steps)
  //   THEN this kills it ...
  //
  if (State(fParticleIsLooping))
  {
    G4double endEnergy = State(fTransportEndKineticEnergy);

    if ((endEnergy < fThreshold_Important_Energy) || (State(fNoLooperTrials)
        >= fThresholdTrials))
    {
      // Kill the looping particle
      //
      // G4cout << "G4ITTransportation will killed the molecule"<< G4endl;
      fParticleChange.ProposeTrackStatus(fStopAndKill);

      // 'Bare' statistics
      fSumEnergyKilled += endEnergy;
      if (endEnergy > fMaxEnergyKilled)
      {
        fMaxEnergyKilled = endEnergy;
      }

#ifdef G4VERBOSE
      if ((fVerboseLevel > 1) || (endEnergy > fThreshold_Warning_Energy))
      {
        G4cout
            << " G4ITTransportation is killing track that is looping or stuck "
            << G4endl<< "   This track has " << track.GetKineticEnergy() / MeV
        << " MeV energy." << G4endl;
        G4cout << "   Number of trials = " << State(fNoLooperTrials)
        << "   No of calls to AlongStepDoIt = " << noCalls
        << G4endl;
      }
#endif
      State(fNoLooperTrials) = 0;
    }
    else
    {
      State(fNoLooperTrials)++;
#ifdef G4VERBOSE
      if ((fVerboseLevel > 2))
      {
        G4cout << "   G4ITTransportation::AlongStepDoIt(): Particle looping -  "
               << "   Number of trials = " << State(fNoLooperTrials)
               << "   No of calls to  = " << noCalls << G4endl;
      }
#endif
    }
  }
  else
  {
    State(fNoLooperTrials)=0;
  }

  // Another (sometimes better way) is to use a user-limit maximum Step size
  // to alleviate this problem ..

  // Introduce smooth curved trajectories to particle-change
  //
  fParticleChange.SetPointerToVectorOfAuxiliaryPoints(
      fFieldPropagator->GimmeTrajectoryVectorAndForgetIt());

#if defined (DEBUG_MEM)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "\t || MEM || End of G4ITTransportation::AlongStepDoIt, diff is: "
      << mem_diff << G4endl;
#endif

  return &fParticleChange;
}

//////////////////////////////////////////////////////////////////////////
//
//  This ensures that the PostStep action is always called,
//  so that it can do the relocation if it is needed.
//

G4double
G4ITTransportation::
PostStepGetPhysicalInteractionLength(const G4Track&, // track
                                     G4double, // previousStepSize
                                     G4ForceCondition* pForceCond)
{
  *pForceCond = Forced;
  return DBL_MAX; // was kInfinity ; but convention now is DBL_MAX
}

/////////////////////////////////////////////////////////////////////////////
//

G4VParticleChange* G4ITTransportation::PostStepDoIt(const G4Track& track,
                                                    const G4Step&)
{
  //    G4cout << "G4ITTransportation::PostStepDoIt" << G4endl;

  PrepareState();
  G4TouchableHandle retCurrentTouchable; // The one to return
  G4bool isLastStep = false;

  // Initialize ParticleChange  (by setting all its members equal
  //                             to corresponding members in G4Track)
  fParticleChange.Initialize(track); // To initialise TouchableChange

  fParticleChange.ProposeTrackStatus(track.GetTrackStatus());

  // If the Step was determined by the volume boundary,
  // logically relocate the particle

  if (State(fGeometryLimitedStep))
  {

    if(fVerboseLevel)
    {
     G4cout << "Step is limited by geometry "
            <<  "track ID : " << track.GetTrackID() << G4endl;
    }

    // fCurrentTouchable will now become the previous touchable,
    // and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)

    if ( State(fCurrentTouchableHandle)->GetVolume() == 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No current touchable found ";
      G4Exception(" G4ITTransportation::PostStepDoIt", "G4ITTransportation001",
                  FatalErrorInArgument, exceptionDescription);
    }

    fLinearNavigator->SetGeometricallyLimitedStep();
    fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle(
        track.GetPosition(), track.GetMomentumDirection(),
        State(fCurrentTouchableHandle), true);
    // Check whether the particle is out of the world volume
    // If so it has exited and must be killed.
    //
    if ( State(fCurrentTouchableHandle)->GetVolume() == 0)
    {
    //  abort();
#ifdef G4VERBOSE
      if (fVerboseLevel > 0)
      {
        G4cout << "Track position : " << track.GetPosition() / nanometer
               << " [nm]" << " Track ID : " << track.GetTrackID() << G4endl;
        G4cout << "G4ITTransportation will killed the track because "
            "State(fCurrentTouchableHandle)->GetVolume() == 0"<< G4endl;
      }
#endif
      fParticleChange.ProposeTrackStatus( fStopAndKill );
    }

    retCurrentTouchable = State(fCurrentTouchableHandle);

//        G4cout << "Current volume : " << track.GetVolume()->GetName()
//               << " Next volume : "
//               << (State(fCurrentTouchableHandle)->GetVolume() ?
//    State(fCurrentTouchableHandle)->GetVolume()->GetName():"OutWorld")
//               << " Position : " << track.GetPosition() / nanometer
//               << " track ID : " << track.GetTrackID()
//               << G4endl;

    fParticleChange.SetTouchableHandle(State(fCurrentTouchableHandle));

    // Update the Step flag which identifies the Last Step in a volume
    isLastStep = fLinearNavigator->ExitedMotherVolume()
               || fLinearNavigator->EnteredDaughterVolume();

#ifdef G4DEBUG_TRANSPORT
    //  Checking first implementation of flagging Last Step in Volume
    G4bool exiting = fLinearNavigator->ExitedMotherVolume();
    G4bool entering = fLinearNavigator->EnteredDaughterVolume();

    if( ! (exiting || entering) )
    {
      G4cout << " Transport> :  Proposed isLastStep= " << isLastStep
      << " Exiting " << fLinearNavigator->ExitedMotherVolume()
      << " Entering " << fLinearNavigator->EnteredDaughterVolume()
      << " Track position : " << track.GetPosition() /nanometer << " [nm]"
      << G4endl;
      G4cout << " Track position : " << track.GetPosition() /nanometer
          << G4endl;
    }
#endif
  }
  else // fGeometryLimitedStep  is false
  {
    // This serves only to move the Navigator's location
    //
//    abort();
    fLinearNavigator->LocateGlobalPointWithinVolume(track.GetPosition());

    // The value of the track's current Touchable is retained.
    // (and it must be correct because we must use it below to
    // overwrite the (unset) one in particle change)
    //  It must be fCurrentTouchable too ??
    //
    fParticleChange.SetTouchableHandle(track.GetTouchableHandle());
    retCurrentTouchable = track.GetTouchableHandle();

    isLastStep = false;
#ifdef G4DEBUG_TRANSPORT
    //  Checking first implementation of flagging Last Step in Volume
    G4cout << " Transport> Proposed isLastStep= " << isLastStep
    << " Geometry did not limit step. Position : "
    << track.GetPosition()/ nanometer << G4endl;
#endif
  } // endif ( fGeometryLimitedStep )

  fParticleChange.ProposeLastStepInVolume(isLastStep);

  const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume();
  const G4Material* pNewMaterial = 0;
  G4VSensitiveDetector* pNewSensitiveDetector = 0;

  if (pNewVol != 0)
  {
    pNewMaterial = pNewVol->GetLogicalVolume()->GetMaterial();
    pNewSensitiveDetector = pNewVol->GetLogicalVolume()->GetSensitiveDetector();
  }

  // ( <const_cast> pNewMaterial ) ;

  fParticleChange.SetMaterialInTouchable((G4Material *) pNewMaterial);
  fParticleChange.SetSensitiveDetectorInTouchable(pNewSensitiveDetector);

  const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
  if (pNewVol != 0)
  {
    pNewMaterialCutsCouple =
        pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();
  }

  if (pNewVol != 0 && pNewMaterialCutsCouple != 0
      && pNewMaterialCutsCouple->GetMaterial() != pNewMaterial)
  {
    // for parametrized volume
    //
    pNewMaterialCutsCouple = G4ProductionCutsTable::GetProductionCutsTable()
        ->GetMaterialCutsCouple(pNewMaterial,
                                pNewMaterialCutsCouple->GetProductionCuts());
  }
  fParticleChange.SetMaterialCutsCoupleInTouchable(pNewMaterialCutsCouple);

  // temporarily until Get/Set Material of ParticleChange,
  // and StepPoint can be made const.
  // Set the touchable in ParticleChange
  // this must always be done because the particle change always
  // uses this value to overwrite the current touchable pointer.
  //
  fParticleChange.SetTouchableHandle(retCurrentTouchable);

  return &fParticleChange;
}

// New method takes over the responsibility to reset the state of
// G4Transportation object at the start of a new track or the resumption of
// a suspended track.

void G4ITTransportation::StartTracking(G4Track* track)
{
  G4VProcess::StartTracking(track);
  if (fInstantiateProcessState)
  {
//        G4VITProcess::fpState = new G4ITTransportationState();
    G4VITProcess::fpState.reset(new G4ITTransportationState());
    // Will set in the same time fTransportationState
  }

  fpSafetyHelper->NewTrackState();
  fpSafetyHelper->SaveTrackState(
      GetIT(track)->GetTrackingInfo()->GetTrackStateManager());

  // The actions here are those that were taken in AlongStepGPIL
  //   when track.GetCurrentStepNumber()==1

  // reset safety value and center
  //
  //    State(fPreviousSafety)    = 0.0 ;
  //    State(fPreviousSftOrigin) = G4ThreeVector(0.,0.,0.) ;

  // reset looping counter -- for motion in field
  //    State(fNoLooperTrials)= 0;
  // Must clear this state .. else it depends on last track's value
  //  --> a better solution would set this from state of suspended track TODO ?
  // Was if( aTrack->GetCurrentStepNumber()==1 ) { .. }

  // ChordFinder reset internal state
  //
  if (DoesGlobalFieldExist())
  {
    fFieldPropagator->ClearPropagatorState();
    // Resets all state of field propagator class (ONLY)
    // including safety values (in case of overlaps and to wipe for first track).

    // G4ChordFinder* chordF= fFieldPropagator->GetChordFinder();
    // if( chordF ) chordF->ResetStepEstimate();
  }

  // Make sure to clear the chord finders of all fields (ie managers)
  static G4ThreadLocal G4FieldManagerStore* fieldMgrStore = 0;
  if (!fieldMgrStore) fieldMgrStore = G4FieldManagerStore::GetInstance();
  fieldMgrStore->ClearAllChordFindersState();

  // Update the current touchable handle  (from the track's)
  //
  PrepareState();
  State(fCurrentTouchableHandle) = track->GetTouchableHandle();

  G4VITProcess::StartTracking(track);
}

#undef State
#undef PrepareState
