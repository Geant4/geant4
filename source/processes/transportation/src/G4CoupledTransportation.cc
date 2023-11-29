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
// ------------------------------------------------------------
//  GEANT 4 class implementation
//
// =======================================================================
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "G4CoupledTransportation.hh"
#include "G4TransportationProcessType.hh"
#include "G4TransportationLogger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ChordFinder.hh"
#include "G4Field.hh"
#include "G4FieldTrack.hh"
#include "G4FieldManagerStore.hh"
#include "G4PathFinder.hh"

#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"

class G4VSensitiveDetector;

G4bool G4CoupledTransportation::fSignifyStepInAnyVolume= false;
// This mode must apply to all threads 

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4CoupledTransportation::G4CoupledTransportation( G4int verbosity )
  : G4Transportation( verbosity, "CoupledTransportation" ),
    fPreviousMassSafety( 0.0 ),
    fPreviousFullSafety( 0.0 ),
    fMassGeometryLimitedStep( false ), 
    fFirstStepInMassVolume( true )
{
  SetProcessSubType(static_cast<G4int>(COUPLED_TRANSPORTATION));
  // SetVerboseLevel is called in the constructor of G4Transportation

  if( verboseLevel > 0 )
  {
    G4cout << " G4CoupledTransportation constructor: ----- " << G4endl;
    G4cout << " Verbose level is " << verboseLevel << G4endl;
    G4cout << " Reports First/Last in " 
           << (fSignifyStepInAnyVolume ? " any " : " mass " )
           << " geometry " << G4endl;
  }
  fPathFinder=  G4PathFinder::GetInstance(); 
}

//////////////////////////////////////////////////////////////////////////

G4CoupledTransportation::~G4CoupledTransportation()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.

G4double G4CoupledTransportation::
AlongStepGetPhysicalInteractionLength( const G4Track&  track,
                                             G4double, //  previousStepSize
                                             G4double  currentMinimumStep,
                                             G4double& proposedSafetyForStart,
                                             G4GPILSelection* selection )
{
  G4double geometryStepLength; 
  G4double startFullSafety= 0.0; // estimated safety for start point (all geometries)
  G4double safetyProposal= -1.0; // local copy of proposal 

  G4ThreeVector  EndUnitMomentum ;
  G4double       lengthAlongCurve = 0.0 ;
 
  fParticleIsLooping = false ;

  // Initial actions moved to  StartTrack()   
  // --------------------------------------
  // Note: in case another process changes touchable handle
  //    it will be necessary to add here (for all steps)   
  // fCurrentTouchableHandle = aTrack->GetTouchableHandle();

  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  //
  *selection = CandidateForSelection ;

  fFirstStepInMassVolume = fNewTrack || fMassGeometryLimitedStep ; 
  fFirstStepInVolume     = fNewTrack || fGeometryLimitedStep ;

#ifdef G4DEBUG_TRANSPORT
  G4cout << "  CoupledTransport::AlongStep GPIL:  "
         << "  1st-step:  any= "  <<fFirstStepInVolume  << " ( geom= "
         << fGeometryLimitedStep << " ) "
         <<           " mass= " << fFirstStepInMassVolume << " ( geom= "
         << fMassGeometryLimitedStep << " ) " 
         << "  newTrack= " << fNewTrack << G4endl;
#endif
  
  // fLastStepInVolume= false;
  fNewTrack = false;

  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
  const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
  G4ThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  G4ThreeVector startPosition          = track.GetPosition() ;
  G4VPhysicalVolume* currentVolume= track.GetVolume(); 

#ifdef G4DEBUG_TRANSPORT
  if( verboseLevel > 1 )
  {
    G4cout << "G4CoupledTransportation::AlongStepGPIL> called in volume " 
           << currentVolume->GetName() << G4endl; 
  }
#endif
  // G4double   theTime        = track.GetGlobalTime() ;

  // The Step Point safety can be limited by other geometries and/or the 
  // assumptions of any process - it's not always the geometrical safety.
  // We calculate the starting point's isotropic safety here.
  //
  G4ThreeVector OriginShift = startPosition - fPreviousSftOrigin ;
  G4double      MagSqShift  = OriginShift.mag2() ;
  startFullSafety= 0.0; 

  // Recall that FullSafety <= MassSafety 
  // Original: if( MagSqShift < sqr(fPreviousMassSafety) ) {
  if( MagSqShift < sqr(fPreviousFullSafety) )
  {
     G4double mag_shift= std::sqrt(MagSqShift); 
     startFullSafety = std::max( (fPreviousFullSafety - mag_shift), 0.0);
       // Need to be consistent between full safety with Mass safety
       // in order reproduce results in simple case
       // --> use same calculation method

     // Only compute full safety if massSafety > 0.  Else it remains 0
     // startFullSafety = fPathFinder->ComputeSafety( startPosition ); 
  }

  // Is the particle charged or has it a magnetic moment?
  //
  G4double particleCharge = pParticle->GetCharge() ;
  G4double magneticMoment = pParticle->GetMagneticMoment() ;
  G4double       restMass = pParticle->GetMass() ; 

  fMassGeometryLimitedStep = false ; //  Set default - alex
  fGeometryLimitedStep     = false;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"

  // Check if the particle has a force, EM or gravitational, exerted on it
  //
  G4FieldManager* fieldMgr= nullptr;
  G4bool          fieldExertsForce = false ;

  const G4Field* ptrField= nullptr;

  fieldMgr = fFieldPropagator->FindAndSetFieldManager( track.GetVolume() );
  G4bool eligibleEM = (particleCharge != 0.0)
                   || ( fUseMagneticMoment && (magneticMoment != 0.0) );
  G4bool eligibleGrav =  fUseGravity && (restMass != 0.0) ;

  if( (fieldMgr!=nullptr) && (eligibleEM||eligibleGrav) )
  {
     // Message the field Manager, to configure it for this track
     //
     fieldMgr->ConfigureForTrack( &track );

     // The above call can transition from a null field-ptr oto a finite field.
     // If the field manager has no field ptr, the field is zero 
     // by definition ( = there is no field ! )
     //
     ptrField= fieldMgr->GetDetectorField();
 
     if( ptrField != nullptr)
     { 
        fieldExertsForce = eligibleEM
              || ( eligibleGrav && ptrField->IsGravityActive() );
     }
  }
  G4double momentumMagnitude = pParticle->GetTotalMomentum() ;

  if( fieldExertsForce )
  {
     auto equationOfMotion= fFieldPropagator->GetCurrentEquationOfMotion();
 
     G4ChargeState chargeState(particleCharge, // Charge can change (dynamic)
                               magneticMoment,
                               pParticleDef->GetPDGSpin() );
     if( equationOfMotion )
     {
        equationOfMotion->SetChargeMomentumMass( chargeState,
                                                 momentumMagnitude,
                                                 restMass );
     }
#ifdef G4DEBUG_TRANSPORT
     else
     {
        G4cerr << " ERROR in G4CoupledTransportation> "
               << "Cannot find valid Equation of motion: "      << G4endl;
               << " Unable to pass Charge, Momentum and Mass "  << G4endl;
     }
#endif     
  }

  G4ThreeVector polarizationVec  = track.GetPolarization() ;
  G4FieldTrack  aFieldTrack = G4FieldTrack(startPosition, 
                                           track.GetGlobalTime(), // Lab.
                                           track.GetMomentumDirection(),
                                           track.GetKineticEnergy(),
                                           restMass,
                                           particleCharge, 
                                           polarizationVec, 
                                           pParticleDef->GetPDGMagneticMoment(),
                                           0.0,  // Length along track
                                           pParticleDef->GetPDGSpin() ) ;
  G4int stepNo= track.GetCurrentStepNumber(); 

  ELimited limitedStep; 
  G4FieldTrack endTrackState('a');  //  Default values

  fMassGeometryLimitedStep = false ;    //  default
  fGeometryLimitedStep     = false;
  if( currentMinimumStep > 0 )
  {
      G4double newMassSafety= 0.0;     //  temp. for recalculation

      // Do the Transport in the field (non recti-linear)
      //
      lengthAlongCurve = fPathFinder->ComputeStep( aFieldTrack,
                                                   currentMinimumStep, 
                                                   G4TransportationManager::kMassNavigatorId,
                                                   stepNo,
                                                   newMassSafety,
                                                   limitedStep,
                                                   endTrackState,
                                                   currentVolume ) ;

      G4double newFullSafety= fPathFinder->GetCurrentSafety();  
        // this was estimated already in step above

      if( limitedStep == kUnique || limitedStep == kSharedTransport )
      {
        fMassGeometryLimitedStep = true ;
      }
      
      fGeometryLimitedStep = (fPathFinder->GetNumberGeometriesLimitingStep() != 0);

#ifdef G4DEBUG_TRANSPORT
      if( fMassGeometryLimitedStep && !fGeometryLimitedStep )
      {
        std::ostringstream message;
        message << " ERROR in determining geometries limiting the step" << G4endl;
        message << "  Limiting:  mass=" << fMassGeometryLimitedStep
                << " any= " << fGeometryLimitedStep << G4endl;
        message << "Incompatible conditions - by which geometries was it limited ?";
        G4Exception("G4CoupledTransportation::AlongStepGetPhysicalInteractionLength()", 
                    "PathFinderConfused", FatalException, message); 
      }
#endif

      geometryStepLength = std::min( lengthAlongCurve, currentMinimumStep); 

      // Momentum:  Magnitude and direction can be changed too now ...
      //
      fMomentumChanged         = true ; 
      fTransportEndMomentumDir = endTrackState.GetMomentumDir() ;

      // Remember last safety origin & value.
      fPreviousSftOrigin  = startPosition ;
      fPreviousMassSafety = newMassSafety ;         
      fPreviousFullSafety = newFullSafety ; 
      // fpSafetyHelper->SetCurrentSafety( newFullSafety, startPosition);

#ifdef G4DEBUG_TRANSPORT
      if( verboseLevel > 1 )
      {
        G4cout << "G4Transport:CompStep> " 
               << " called the pathfinder for a new step at " << startPosition
               << " and obtained step = " << lengthAlongCurve << G4endl;
        G4cout << "  New safety (preStep) = " << newMassSafety << G4endl;
      }
#endif

      // Store as best estimate value
      startFullSafety    = newFullSafety ; 

      // Get the End-Position and End-Momentum (Dir-ection)
      fTransportEndPosition = endTrackState.GetPosition() ;
      fTransportEndKineticEnergy  = endTrackState.GetKineticEnergy() ; 
  }
  else
  {
      geometryStepLength   = lengthAlongCurve= 0.0 ;
      fMomentumChanged         = false ; 
      // fMassGeometryLimitedStep = false ;   //  --- ???
      // fGeometryLimitedStep = true;
      fTransportEndMomentumDir = track.GetMomentumDirection();
      fTransportEndKineticEnergy  = track.GetKineticEnergy();

      fTransportEndPosition = startPosition;

      endTrackState= aFieldTrack;  // Ensures that time is updated
  }
  // G4FieldTrack aTrackState(endTrackState);  

  if( !fieldExertsForce ) 
  { 
      fParticleIsLooping         = false ; 
      fMomentumChanged           = false ; 
      fEndGlobalTimeComputed     = false ; 
  } 
  else 
  { 
      fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;
  
#ifdef G4DEBUG_TRANSPORT
      if( verboseLevel > 1 )
      {
        G4cout << " G4CT::CS End Position = "
               << fTransportEndPosition << G4endl; 
        G4cout << " G4CT::CS End Direction = "
               << fTransportEndMomentumDir << G4endl; 
      }
#endif
      if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() )
      {
          // If the field can change energy, then the time must be integrated
          //    - so this should have been updated
          //
          fCandidateEndGlobalTime   = endTrackState.GetLabTimeOfFlight();
          fEndGlobalTimeComputed    = true;
  
          // was ( fCandidateEndGlobalTime != track.GetGlobalTime() );
          // a cleaner way is to have FieldTrack knowing whether time
          // is updated
      }
      else
      {
          // The energy should be unchanged by field transport,
          //    - so the time changed will be calculated elsewhere
          //
          fEndGlobalTimeComputed = false;
  
#ifdef G4VERBOSE
          // Check that the integration preserved the energy 
          //     -  and if not correct this!
          G4double  startEnergy= track.GetKineticEnergy();
          G4double  endEnergy= fTransportEndKineticEnergy; 
      
          G4double absEdiff = std::fabs(startEnergy- endEnergy);
          if( (verboseLevel > 1) && ( absEdiff > perThousand * endEnergy) )
          {
            ReportInexactEnergy(startEnergy, endEnergy); 
          }  // end of if (verboseLevel)
#endif
          // Correct the energy for fields that conserve it
          //  This - hides the integration error
          //       - but gives a better physical answer
          fTransportEndKineticEnergy= track.GetKineticEnergy();
      }
  }

  fEndPointDistance   = (fTransportEndPosition - startPosition).mag() ;
  fTransportEndSpin = endTrackState.GetSpin();

  // Calculate the safety

  safetyProposal= startFullSafety;   // used to be startMassSafety
    // Changed to accomodate processes that cannot update the safety

  // Update safety for the end-point, if becomes negative at the end-point.

  if(   (startFullSafety < fEndPointDistance )
        && ( particleCharge != 0.0 ) )  // Only needed to prepare for MSC
   //   && !fGeometryLimitedStep ) // To-Try: No safety update if at boundary
  {
      G4double endFullSafety =
        fPathFinder->ComputeSafety( fTransportEndPosition); 
        // Expected mission -- only mass geometry's safety
        //   fLinearNavigator->ComputeSafety( fTransportEndPosition) ;
        // Yet discrete processes only have poststep -- and this cannot 
        //   currently revise the safety  
        //   ==> so we use the all-geometry safety as a precaution

      fpSafetyHelper->SetCurrentSafety( endFullSafety, fTransportEndPosition);
        // Pushing safety to Helper avoids recalculation at this point

      G4ThreeVector centerPt= G4ThreeVector(0.0, 0.0, 0.0);  // Used for return value
      G4double endMassSafety= fPathFinder->ObtainSafety( G4TransportationManager::kMassNavigatorId, centerPt);
        //  Retrieves the mass value from PathFinder (it calculated it)

      fPreviousMassSafety = endMassSafety ; 
      fPreviousFullSafety = endFullSafety; 
      fPreviousSftOrigin = fTransportEndPosition ;

      // The convention (Stepping Manager's) is safety from the start point
      //
      safetyProposal = endFullSafety + fEndPointDistance;
          //  --> was endMassSafety
      // Changed to accomodate processes that cannot update the safety

#ifdef G4DEBUG_TRANSPORT 
      G4int prec= G4cout.precision(12) ;
      G4cout << "***CoupledTransportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Revised Safety at endpoint "  << fTransportEndPosition
             << "   give safety values: Mass= " << endMassSafety 
             << "  All= " << endFullSafety << G4endl ; 
      G4cout << "  Adding endpoint distance " << fEndPointDistance
             << "   to obtain pseudo-safety= " << safetyProposal << G4endl ; 
      G4cout.precision(prec); 
  }  
  else
  {
      G4int prec= G4cout.precision(12) ;
      G4cout << "***CoupledTransportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Quick Safety estimate at endpoint "
             << fTransportEndPosition
             << "   gives safety endpoint value = "
             << startFullSafety - fEndPointDistance
             << "  using start-point value " << startFullSafety 
             << "  and endpointDistance " << fEndPointDistance << G4endl;
      G4cout.precision(prec); 
#endif
  }          

  proposedSafetyForStart= safetyProposal; 
  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

/////////////////////////////////////////////////////////////////////////////

void G4CoupledTransportation::
ReportMove( G4ThreeVector OldVector, G4ThreeVector NewVector,
            const G4String& Quantity )
{
    G4ThreeVector moveVec = ( NewVector - OldVector );

    G4cerr << G4endl
           << "**************************************************************"
           << G4endl;
    G4cerr << "Endpoint has moved between value expected from TransportEndPosition "
           << " and value from Track in PostStepDoIt. " << G4endl
           << "Change of " << Quantity << " is " << moveVec.mag() / mm
           << " mm long, "
           << " and its vector is " << (1.0/mm) * moveVec << " mm " << G4endl
           << "Endpoint of ComputeStep was " << OldVector
           << " and current position to locate is " << NewVector << G4endl;
}

/////////////////////////////////////////////////////////////////////////////

G4VParticleChange* G4CoupledTransportation::PostStepDoIt( const G4Track& track,
                                                          const G4Step& )
{
  G4TouchableHandle retCurrentTouchable ;   // The one to return

  // Initialize ParticleChange  (by setting all its members equal
  //                             to corresponding members in G4Track)
  // fParticleChange.Initialize(track) ;  // To initialise TouchableChange

  fParticleChange.ProposeTrackStatus(track.GetTrackStatus()) ;

  if( fSignifyStepInAnyVolume )
  {
     fParticleChange.ProposeFirstStepInVolume( fFirstStepInVolume );
  }
  else
  {
     fParticleChange.ProposeFirstStepInVolume( fFirstStepInMassVolume );
  }
  
  // Check that the end position and direction are preserved 
  // since call to AlongStepDoIt

#ifdef G4DEBUG_TRANSPORT
  if( ( verboseLevel > 0 )
     && ((fTransportEndPosition - track.GetPosition()).mag2() >= 1.0e-16) )
  {
     ReportMove( track.GetPosition(), fTransportEndPosition,
                 "End of Step Position" ); 
     G4cerr << " Problem in G4CoupledTransportation::PostStepDoIt " << G4endl; 
  }

  // If the Step was determined by the volume boundary, relocate the particle
  // The pathFinder will know that the geometry limited the step (!?)

  if( verboseLevel > 0 )
  {
     G4cout << " Calling PathFinder::Locate() from " 
            << " G4CoupledTransportation::PostStepDoIt " << G4endl;
     G4cout << "   fGeometryLimitedStep is " << fGeometryLimitedStep << G4endl;
  }
#endif

  if(fGeometryLimitedStep)
  {  
    fPathFinder->Locate( track.GetPosition(), 
                         track.GetMomentumDirection(),
                         true); 

    // fCurrentTouchable will now become the previous touchable, 
    // and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)

    fCurrentTouchableHandle= 
      fPathFinder->CreateTouchableHandle( G4TransportationManager::kMassNavigatorId );

#ifdef G4DEBUG_TRANSPORT
    if( verboseLevel > 1 )
    {
       G4VPhysicalVolume* vol= fCurrentTouchableHandle->GetVolume(); 
       G4cout << "CHECK !!!!!!!!!!! fCurrentTouchableHandle->GetVolume() = "
              << vol;
       if( vol ) { G4cout << "Name=" << vol->GetName(); }
       G4cout << G4endl;
    }
#endif

    // Check whether the particle is out of the world volume 
    // If so it has exited and must be killed.
    //
    if( fCurrentTouchableHandle->GetVolume() == 0 )
    {
       fParticleChange.ProposeTrackStatus( fStopAndKill ) ;
    }
    retCurrentTouchable = fCurrentTouchableHandle ;
    // fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;
  }
  else  // fGeometryLimitedStep  is false
  { 
#ifdef G4DEBUG_TRANSPORT
    if( verboseLevel > 1 )
    {
      G4cout << "G4CoupledTransportation::PostStepDoIt -- "
             << " fGeometryLimitedStep  = " << fGeometryLimitedStep
             << " must be false " << G4endl;
    }
#endif
    // This serves only to move each of the Navigator's location
    //
    // fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

    fPathFinder->ReLocate( track.GetPosition() );
                           // track.GetMomentumDirection() ); 

    // Keep the value of the track's current Touchable is retained,
    //  and use it to overwrite the (unset) one in particle change.
    // Expect this must be fCurrentTouchable too
    //   - could it be different, eg at the start of a step ?
    //
    retCurrentTouchable = track.GetTouchableHandle() ;
    // fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
  }  // endif ( fGeometryLimitedStep )

#ifdef G4DEBUG_NAVIGATION  
  G4cout << "  CoupledTransport::AlongStep GPIL:  "
         << " last-step:  any= " << fGeometryLimitedStep << " . ..... x . "
         << " mass= " << fMassGeometryLimitedStep << G4endl;
#endif
  
  if( fSignifyStepInAnyVolume )
    fParticleChange.ProposeLastStepInVolume(fGeometryLimitedStep);
  else
     fParticleChange.ProposeLastStepInVolume(fMassGeometryLimitedStep);
  
  SetTouchableInformation(retCurrentTouchable);

  return &fParticleChange ;
}

/////////////////////////////////////////////////////////////////////////////
// New method takes over the responsibility to reset the state of 
// G4CoupledTransportation object:
//      - at the start of a new track,  and
//      - on the resumption of a suspended track. 
//
void 
G4CoupledTransportation::StartTracking(G4Track* aTrack)
{
  G4Transportation::StartTracking(aTrack);
  
  G4ThreeVector position = aTrack->GetPosition(); 
  G4ThreeVector direction = aTrack->GetMomentumDirection();

  fPathFinder->PrepareNewTrack( position, direction); 
  // This implies a call to fPathFinder->Locate( position, direction ); 

  // reset safety value and center
  //
  fPreviousMassSafety  = 0.0 ; 
  fPreviousFullSafety  = 0.0 ; 
  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
}

/////////////////////////////////////////////////////////////////////////////

void 
G4CoupledTransportation::EndTracking()
{
  G4TransportationManager::GetTransportationManager()->InactivateAll();
  fPathFinder->EndTrack(); 
    // Resets TransportationManager to use ordinary Navigator
}

/////////////////////////////////////////////////////////////////////////////

void
G4CoupledTransportation::
ReportInexactEnergy(G4double startEnergy, G4double endEnergy)
{
  static G4ThreadLocal G4int no_warnings= 0, warnModulo=1,
                             moduloFactor= 10, no_large_ediff= 0; 

  if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
  {
    no_large_ediff ++;
    if( (no_large_ediff% warnModulo) == 0 )
    {
      no_warnings++;
      std::ostringstream message;
      message << "Energy change in Step is above 1^-3 relative value. "
              << G4endl
              << "   Relative change in 'tracking' step = " 
              << std::setw(15) << (endEnergy-startEnergy)/startEnergy
              << G4endl
              << "   Starting E= " << std::setw(12) << startEnergy / MeV
              << " MeV " << G4endl
              << "   Ending   E= " << std::setw(12) << endEnergy   / MeV
              << " MeV " << G4endl
              << "Energy has been corrected -- however, review"
              << " field propagation parameters for accuracy." << G4endl;
      if ( (verboseLevel > 2 ) || (no_warnings<4)
        || (no_large_ediff == warnModulo * moduloFactor) )
      {
        message << "These include EpsilonStepMax(/Min) in G4FieldManager,"
                << G4endl
                << "which determine fractional error per step for integrated quantities."
                << G4endl
                << "Note also the influence of the permitted number of integration steps."
                << G4endl;
      }
      message << "Bad 'endpoint'. Energy change detected and corrected."
              << G4endl
              << "Has occurred already " << no_large_ediff << " times.";
      G4Exception("G4CoupledTransportation::AlongStepGetPIL()", 
                  "EnergyChange", JustWarning, message);
      if( no_large_ediff == warnModulo * moduloFactor )
      {
        warnModulo *= moduloFactor;
      }
    }
  }
}
