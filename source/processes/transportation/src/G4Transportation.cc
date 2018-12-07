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
//  GEANT 4  include file implementation
//
// ------------------------------------------------------------
//
// This class is a process responsible for the transportation of 
// a particle, ie the geometrical propagation that encounters the 
// geometrical sub-volumes of the detectors.
//
// It is also tasked with the key role of proposing the "isotropic safety",
//   which will be used to update the post-step point's safety.
//
// =======================================================================
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "G4Transportation.hh"
#include "G4TransportationProcessType.hh"
#include "G4TransportationLogger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"

#include "G4FieldManagerStore.hh"
#include "G4CoupledTransportation.hh"

class G4VSensitiveDetector;

G4bool G4Transportation::fUseMagneticMoment=false;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4Transportation::G4Transportation( G4int verbosity )
  : G4VProcess( G4String("Transportation"), fTransportation ),
    // fTransportEndPosition( 0.0, 0.0, 0.0 ),
    // fTransportEndMomentumDir( 0.0, 0.0, 0.0 ),
    // fTransportEndKineticEnergy( 0.0 ),
    // fTransportEndSpin( 0.0, 0.0, 0.0 ),
    // fMomentumChanged(true),
    // fEndGlobalTimeComputed(false), 
    // fCandidateEndGlobalTime(0.0),
    // fParticleIsLooping( false ),
    // fNewTrack( true ),
    // fFirstStepInVolume( true ),
    // fLastStepInVolume( false ), 
    // fGeometryLimitedStep(true),
    fFieldExertedForce( false ),
    fPreviousSftOrigin( 0.,0.,0. ),
    fPreviousSafety( 0.0 ),
    fEndPointDistance( -1.0 ), 
    fShortStepOptimisation( false ) // Old default: true (=fast short steps)
{
  SetProcessSubType(static_cast<G4int>(TRANSPORTATION));
  pParticleChange= &fParticleChange;   // Required to conform to G4VProcess 
  SetVerboseLevel(verbosity);

  G4TransportationManager* transportMgr ; 

  transportMgr = G4TransportationManager::GetTransportationManager() ; 

  fLinearNavigator = transportMgr->GetNavigatorForTracking() ; 

  fFieldPropagator = transportMgr->GetPropagatorInField() ;

  fpSafetyHelper =   transportMgr->GetSafetyHelper();  // New 

  fpLogger = new G4TransportationLogger("G4Transportation", verbosity);

  SetHighLooperThresholds();
  // Use the old defaults: Warning = 100 MeV, Important = 250 MeV, No Trials = 10;
  
  PushThresholdsToLogger();
  // Should be done by Set methods in SetHighLooperThresholds -- making sure
  
  // Cannot determine whether a field exists here, as it would 
  //  depend on the relative order of creating the detector's 
  //  field and this process. That order is not guaranted.
  // Instead later the method DoesGlobalFieldExist() is called

  static G4ThreadLocal G4TouchableHandle* pNullTouchableHandle = 0;
  if ( !pNullTouchableHandle)
  {
    pNullTouchableHandle = new G4TouchableHandle;
  }
  fCurrentTouchableHandle = *pNullTouchableHandle;
    // Points to (G4VTouchable*) 0

#ifdef G4VERBOSE
  if( verboseLevel > 0) 
  { 
     G4cout << " G4Transportation constructor> set fShortStepOptimisation to "; 
     if ( fShortStepOptimisation )  { G4cout << "true"  << G4endl; }
     else                           { G4cout << "false" << G4endl; }
  } 
#endif
}

//////////////////////////////////////////////////////////////////////////

G4Transportation::~G4Transportation()
{
  if( fSumEnergyKilled > 0.0 )
  {
     PrintStatistics( G4cout );     
  }
  delete fpLogger;
}

//////////////////////////////////////////////////////////////////////////

void
G4Transportation::PrintStatistics( std::ostream& outStr) const
{
   outStr << " G4Transportation: Statistics for looping particles " << G4endl;   
   if( fSumEnergyKilled > 0.0 || fNumLoopersKilled > 0 )
   {
      outStr << "   Sum of energy of looping tracks killed: "
             <<  fSumEnergyKilled / CLHEP::MeV << " MeV "          
             << " from " << fNumLoopersKilled << "  tracks "  << G4endl
             <<  "  Sum of energy of non-electrons        : "
             << fSumEnergyKilled_NonElectron / CLHEP::MeV << " MeV "
             << "  from " << fNumLoopersKilled_NonElectron << " tracks "
             << G4endl;
      outStr << "   Max energy of  *any type*  looper killed: " << fMaxEnergyKilled
             << "    its PDG was " << fMaxEnergyKilledPDG  << G4endl;
      if( fMaxEnergyKilled_NonElectron > 0.0 )
      {
         outStr << "   Max energy of non-electron looper killed: " 
                << fMaxEnergyKilled_NonElectron
                << "    its PDG was " << fMaxEnergyKilled_NonElecPDG << G4endl;
      }
      if( fMaxEnergySaved > 0.0 )
      {
         outStr << "   Max energy of loopers 'saved':  " << fMaxEnergySaved << G4endl;
         outStr << "   Sum of energy of loopers 'saved': "
                <<  fSumEnergySaved << G4endl;
         outStr << "   Sum of energy of unstable loopers 'saved': "
                << fSumEnergyUnstableSaved << G4endl;
      }
   }
   else
   {
      outStr << " No looping tracks found or killed. " << G4endl;
   }
}

//////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.

G4double G4Transportation::
AlongStepGetPhysicalInteractionLength( const G4Track&  track,
                                             G4double, //  previousStepSize
                                             G4double  currentMinimumStep,
                                             G4double& currentSafety,
                                             G4GPILSelection* selection )
{
  G4double geometryStepLength= -1.0, newSafety= -1.0; 
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

  fFirstStepInVolume= fNewTrack || fLastStepInVolume;
  fLastStepInVolume= false;
  fNewTrack = false;

  fParticleChange.ProposeFirstStepInVolume(fFirstStepInVolume);
  
  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle*    pParticle    = track.GetDynamicParticle() ;
  const G4ParticleDefinition* pParticleDef = pParticle->GetDefinition() ;
  G4ThreeVector startMomentumDir = pParticle->GetMomentumDirection() ;
  G4ThreeVector startPosition    = track.GetPosition() ;

  // The Step Point safety can be limited by other geometries and/or the 
  // assumptions of any process - it's not always the geometrical safety.
  // We calculate the starting point's isotropic safety here.
  //
  G4ThreeVector OriginShift = startPosition - fPreviousSftOrigin ;
  G4double      MagSqShift  = OriginShift.mag2() ;
  if( MagSqShift >= sqr(fPreviousSafety) )
  {
     currentSafety = 0.0 ;
  }
  else
  {
     currentSafety = fPreviousSafety - std::sqrt(MagSqShift) ;
  }

  // Is the particle charged or has it a magnetic moment?
  //
  G4double particleCharge = pParticle->GetCharge() ; 
  G4double magneticMoment = pParticle->GetMagneticMoment() ;
  G4double       restMass = pParticle->GetMass() ;

  fGeometryLimitedStep = false ;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"

  // Check if the particle has a force, EM or gravitational, exerted on it
  //
  G4FieldManager* fieldMgr=0;
  G4bool          fieldExertsForce = false ;

  G4bool gravityOn = false;
  G4bool fieldExists = false;  // Field is not 0 (null pointer)

  fieldMgr = fFieldPropagator->FindAndSetFieldManager( track.GetVolume() );
  if( fieldMgr != 0 )
  {
     // Message the field Manager, to configure it for this track
     //
     fieldMgr->ConfigureForTrack( &track );

     // Is here to allow a transition from no-field pointer 
     // to finite field (non-zero pointer).
     // If the field manager has no field ptr, the field is zero 
     // by definition ( = there is no field ! )
     //
     const G4Field* ptrField= fieldMgr->GetDetectorField();
     fieldExists = (ptrField!=0) ;
     if( fieldExists ) 
     {
        gravityOn= ptrField->IsGravityActive();

        if( (particleCharge != 0.0) 
            || (fUseMagneticMoment && (magneticMoment != 0.0) )
            || (gravityOn          && (restMass != 0.0) ) )
        {
          fieldExertsForce = fieldExists; 
        }
     }
  }
  fFieldExertedForce = fieldExertsForce; 
  
  if( !fieldExertsForce ) 
  {
     G4double linearStepLength ;
     if( fShortStepOptimisation && (currentMinimumStep <= currentSafety) )
     {
       // The Step is guaranteed to be taken
       //
       geometryStepLength   = currentMinimumStep ;
       fGeometryLimitedStep = false ;
     }
     else
     {
       //  Find whether the straight path intersects a volume
       //
       linearStepLength = fLinearNavigator->ComputeStep( startPosition, 
                                                         startMomentumDir,
                                                         currentMinimumStep, 
                                                         newSafety) ;
       // Remember last safety origin & value.
       //
       fPreviousSftOrigin = startPosition ;
       fPreviousSafety    = newSafety ; 
       fpSafetyHelper->SetCurrentSafety( newSafety, startPosition);

       currentSafety = newSafety ;
          
       fGeometryLimitedStep= (linearStepLength <= currentMinimumStep); 
       if( fGeometryLimitedStep )
       {
         // The geometry limits the Step size (an intersection was found.)
         geometryStepLength   = linearStepLength ;
       } 
       else
       {
         // The full Step is taken.
         geometryStepLength   = currentMinimumStep ;
       }
     }
     fEndPointDistance = geometryStepLength ;

     // Calculate final position
     //
     fTransportEndPosition = startPosition+geometryStepLength*startMomentumDir;

     // Momentum direction, energy and polarisation are unchanged by transport
     //
     fTransportEndMomentumDir   = startMomentumDir ; 
     fTransportEndKineticEnergy = track.GetKineticEnergy() ;
     fTransportEndSpin          = track.GetPolarization();
     fParticleIsLooping         = false ;
     fMomentumChanged           = false ; 
     fEndGlobalTimeComputed     = false ;
  }
  else   //  A field exerts force
  {
     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
     G4ThreeVector  EndUnitMomentum ;
     G4double       lengthAlongCurve ;

     // The charge can change (dynamic)
     //
     G4ChargeState chargeState(particleCharge,
                               magneticMoment,
                               pParticleDef->GetPDGSpin() );
     
     auto equationOfMotion = fFieldPropagator->GetCurrentEquationOfMotion();

     equationOfMotion->SetChargeMomentumMass( chargeState,
                                              momentumMagnitude,
                                              restMass);
      
     G4FieldTrack aFieldTrack = G4FieldTrack( startPosition, 
                                           track.GetGlobalTime(), // Lab.
                                           track.GetMomentumDirection(),
                                           track.GetKineticEnergy(),
                                           restMass,
                                           particleCharge, 
                                           track.GetPolarization(),
                                           pParticleDef->GetPDGMagneticMoment(),
                                           0.0,  // Length along track
                                           pParticleDef->GetPDGSpin() );

     if( currentMinimumStep > 0 ) 
     {
        // Do the Transport in the field (non recti-linear)
        //
        lengthAlongCurve = fFieldPropagator->ComputeStep( aFieldTrack,
                                                          currentMinimumStep, 
                                                          currentSafety,
                                                          track.GetVolume() );

        fGeometryLimitedStep= fFieldPropagator->IsLastStepInVolume();
        //
        // It is possible that step was reduced in PropagatorInField due to
        // previous zero steps. To cope with case that reduced step is taken
        // in full, we must rely on PiF to obtain this value

        geometryStepLength = std::min( lengthAlongCurve, currentMinimumStep );
        
        // Remember last safety origin & value.
        //
        fPreviousSftOrigin = startPosition ;
        fPreviousSafety    = currentSafety ;         
        fpSafetyHelper->SetCurrentSafety( currentSafety, startPosition);
     }
     else
     {
        geometryStepLength   = lengthAlongCurve= 0.0 ;
        fGeometryLimitedStep = false ;
     }
      
     // Get the End-Position and End-Momentum (Dir-ection)
     //
     fTransportEndPosition = aFieldTrack.GetPosition() ;

     // Momentum:  Magnitude and direction can be changed too now ...
     //
     fMomentumChanged         = true ; 
     fTransportEndMomentumDir = aFieldTrack.GetMomentumDir() ;

     fTransportEndKineticEnergy  = aFieldTrack.GetKineticEnergy() ; 

     if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() )
     {
        // If the field can change energy, then the time must be integrated
        //    - so this should have been updated
        //
        fCandidateEndGlobalTime   = aFieldTrack.GetLabTimeOfFlight();
        fEndGlobalTimeComputed    = true;

        // was ( fCandidateEndGlobalTime != track.GetGlobalTime() );
        // a cleaner way is to have FieldTrack knowing whether time is updated.
     }
     else
     {
        // The energy should be unchanged by field transport,
        //    - so the time changed will be calculated elsewhere
        //
        fEndGlobalTimeComputed = false;

        // Check that the integration preserved the energy 
        //     -  and if not correct this!
        G4double  startEnergy= track.GetKineticEnergy();
        G4double  endEnergy= fTransportEndKineticEnergy; 

        static G4ThreadLocal G4int no_inexact_steps=0, no_large_ediff;
        G4double absEdiff = std::fabs(startEnergy- endEnergy);
        if( absEdiff > perMillion * endEnergy )
        {
          no_inexact_steps++;
          // Possible statistics keeping here ...
        }
        if( verboseLevel > 1 )
        {
          if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
          {
            static G4ThreadLocal G4int no_warnings= 0, warnModulo=1,
                                       moduloFactor= 10; 
            no_large_ediff ++;
            if( (no_large_ediff% warnModulo) == 0 )
            {
              no_warnings++;
              std::ostringstream message;
              message << "Energy change in Step is above 1^-3 relative value. "
                << G4endl
                << "     Relative change in 'tracking' step = " 
                << std::setw(15) << (endEnergy-startEnergy)/startEnergy
                << G4endl
                << "     Starting E= " << std::setw(12) << startEnergy / MeV
                << " MeV " << G4endl
                << "     Ending   E= " << std::setw(12) << endEnergy   / MeV
                << " MeV " << G4endl       
                << "Energy has been corrected -- however, review"
                << " field propagation parameters for accuracy." << G4endl;
              if ( (verboseLevel > 2 ) || (no_warnings<4)
                || (no_large_ediff == warnModulo * moduloFactor) )
              {
                message << "These include EpsilonStepMax(/Min) in G4FieldManager " << G4endl
                        << "which determine fractional error per step for integrated quantities. " << G4endl
                        << "Note also the influence of the permitted number of integration steps."
                        << G4endl;
              }
              message << "Bad 'endpoint'. Energy change detected and corrected."
                      << G4endl
                      << "Has occurred already " << no_large_ediff << " times.";
              G4Exception("G4Transportation::AlongStepGetPIL()", 
                          "EnergyChange", JustWarning, message); 
              if( no_large_ediff == warnModulo * moduloFactor )
              {
                 warnModulo *= moduloFactor;
              }
            }
          }
        }  // end of if (verboseLevel)

        // Correct the energy for fields that conserve it
        //  This - hides the integration error
        //       - but gives a better physical answer
        //
        fTransportEndKineticEnergy= track.GetKineticEnergy(); 
     }

     fTransportEndSpin = aFieldTrack.GetSpin();
     fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;
     fEndPointDistance   = (fTransportEndPosition - startPosition).mag() ;
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if( currentMinimumStep == 0.0 ) 
  {
      if( currentSafety == 0.0 )  { fGeometryLimitedStep = true; }
  }

  // Update the safety starting from the end-point,
  // if it will become negative at the end-point.
  //
  if( currentSafety < fEndPointDistance ) 
  {
      if( particleCharge != 0.0 ) 
      {
         G4double endSafety =
               fLinearNavigator->ComputeSafety( fTransportEndPosition) ;
         currentSafety      = endSafety ;
         fPreviousSftOrigin = fTransportEndPosition ;
         fPreviousSafety    = currentSafety ; 
         fpSafetyHelper->SetCurrentSafety(currentSafety, fTransportEndPosition);

         // Because the Stepping Manager assumes it is from the start point, 
         //  add the StepLength
         //
         currentSafety     += fEndPointDistance ;

#ifdef G4DEBUG_TRANSPORT 
         G4cout.precision(12) ;
         G4cout << "***G4Transportation::AlongStepGPIL ** " << G4endl  ;
         G4cout << "  Called Navigator->ComputeSafety at " << fTransportEndPosition
                << "    and it returned safety=  " << endSafety << G4endl ; 
         G4cout << "  Adding endpoint distance   " << fEndPointDistance  
                << "    to obtain pseudo-safety= " << currentSafety << G4endl ; 
      }
      else
      {
         G4cout << "***G4Transportation::AlongStepGPIL ** " << G4endl  ;
         G4cout << "  Avoiding call to ComputeSafety : " << G4endl;
         G4cout << "    charge     = " << particleCharge     << G4endl;
         G4cout << "    mag moment = " << magneticMoment     << G4endl;
#endif
      }
  }            

  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)

G4VParticleChange* G4Transportation::AlongStepDoIt( const G4Track& track,
                                                    const G4Step&  stepData )
{
  static G4ThreadLocal G4long noCallsASDI=0;
  const char *methodName= "AlongStepDoIt";
  noCallsASDI++;

  fParticleChange.Initialize(track) ;

  //  Code for specific process 
  //
  fParticleChange.ProposePosition(fTransportEndPosition) ;
  fParticleChange.ProposeMomentumDirection(fTransportEndMomentumDir) ;
  fParticleChange.ProposeEnergy(fTransportEndKineticEnergy) ;
  fParticleChange.SetMomentumChanged(fMomentumChanged) ;

  fParticleChange.ProposePolarization(fTransportEndSpin);

  G4double deltaTime = 0.0 ;

  // Calculate  Lab Time of Flight (ONLY if field Equations used it!)
  // G4double endTime   = fCandidateEndGlobalTime;
  // G4double delta_time = endTime - startTime;

  G4double startTime = track.GetGlobalTime() ;
  
  if (!fEndGlobalTimeComputed)
  {
     // The time was not integrated .. make the best estimate possible
     //
     G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity();
     G4double stepLength      = track.GetStepLength();

     deltaTime= 0.0;  // in case initialVelocity = 0 
     if ( initialVelocity > 0.0 )  { deltaTime = stepLength/initialVelocity; }

     fCandidateEndGlobalTime   = startTime + deltaTime ;
     fParticleChange.ProposeLocalTime(  track.GetLocalTime() + deltaTime) ;
  }
  else
  {
     deltaTime = fCandidateEndGlobalTime - startTime ;
     fParticleChange.ProposeGlobalTime( fCandidateEndGlobalTime ) ;
  }


  // Now Correct by Lorentz factor to get delta "proper" Time
  
  G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
  G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;

  fParticleChange.ProposeProperTime(track.GetProperTime() + deltaProperTime) ;
  //fParticleChange.ProposeTrueStepLength( track.GetStepLength() ) ;

  // If the particle is caught looping or is stuck (in very difficult
  // boundaries) in a magnetic field (doing many steps) THEN can kill it ...
  //
  if ( fParticleIsLooping )
  {
      G4double endEnergy= fTransportEndKineticEnergy;
      fNoLooperTrials ++; 
      auto particleType= track.GetDynamicParticle()->GetParticleDefinition();
      
      G4bool stable = particleType->GetPDGStable();
      G4bool candidateForEnd = (endEnergy < fThreshold_Important_Energy) 
                            || (fNoLooperTrials >= fThresholdTrials) ;
      G4bool unstableAndKillable = !stable && ( fAbandonUnstableTrials != 0);
      G4bool unstableForEnd = (endEnergy < fThreshold_Important_Energy)
                              && (fNoLooperTrials >= fAbandonUnstableTrials) ;
      if( (candidateForEnd && stable) || (unstableAndKillable && unstableForEnd) )
      {
        // Kill the looping particle 
        //
        fParticleChange.ProposeTrackStatus( fStopAndKill )  ;
        G4int particlePDG= particleType->GetPDGEncoding();
        const G4int electronPDG= 11; // G4Electron::G4Electron()->GetPDGEncoding();

        // Simple statistics
        fSumEnergyKilled += endEnergy;
        fSumEnerSqKilled = endEnergy * endEnergy;
        fNumLoopersKilled++;
        
        if( endEnergy > fMaxEnergyKilled ) {
           fMaxEnergyKilled = endEnergy;
           fMaxEnergyKilledPDG = particlePDG; 
        }
        if(  particleType->GetPDGEncoding() != electronPDG )
        {
           fSumEnergyKilled_NonElectron += endEnergy;
           fSumEnerSqKilled_NonElectron += endEnergy * endEnergy;
           fNumLoopersKilled_NonElectron++;
           
           if( endEnergy > fMaxEnergyKilled_NonElectron )
           {
              fMaxEnergyKilled_NonElectron = endEnergy;
              fMaxEnergyKilled_NonElecPDG =  particlePDG;
           }
        }

        if( endEnergy > fThreshold_Warning_Energy )
        {
          fpLogger->ReportLoopingTrack( track, stepData, fNoLooperTrials,
                                        noCallsASDI, methodName );
        }
        fNoLooperTrials=0; 
      }
      else
      {
        fMaxEnergySaved = std::max( endEnergy, fMaxEnergySaved);
        if( fNoLooperTrials == 1 ) {
          fSumEnergySaved += endEnergy;
          if ( !stable )
             fSumEnergyUnstableSaved += endEnergy;
        }
#ifdef G4VERBOSE
        if( verboseLevel > 2 )
        {
          G4cout << "   " << methodName  
                 << " Particle is looping but is saved ..."  << G4endl             
                 << "   Number of trials = " << fNoLooperTrials << G4endl
                 << "   No of calls to  = " << noCallsASDI << G4endl;
        }
#endif
      }
  }
  else
  {
      fNoLooperTrials=0; 
  }

  // Another (sometimes better way) is to use a user-limit maximum Step size
  // to alleviate this problem .. 

  // Introduce smooth curved trajectories to particle-change
  //
  fParticleChange.SetPointerToVectorOfAuxiliaryPoints
    (fFieldPropagator->GimmeTrajectoryVectorAndForgetIt() );

  return &fParticleChange ;
}

//////////////////////////////////////////////////////////////////////////
//
//  This ensures that the PostStep action is always called,
//  so that it can do the relocation if it is needed.
// 

G4double G4Transportation::
PostStepGetPhysicalInteractionLength( const G4Track&,
                                            G4double, // previousStepSize
                                            G4ForceCondition* pForceCond )
{
  fFieldExertedForce = false; // Not known
  *pForceCond = Forced ; 
  return DBL_MAX ;  // was kInfinity ; but convention now is DBL_MAX
}

/////////////////////////////////////////////////////////////////////////////
//

G4VParticleChange* G4Transportation::PostStepDoIt( const G4Track& track,
                                                   const G4Step& )
{
   G4TouchableHandle retCurrentTouchable ;   // The one to return
   G4bool isLastStep= false; 

  // Initialize ParticleChange  (by setting all its members equal
  //                             to corresponding members in G4Track)
  // fParticleChange.Initialize(track) ;  // To initialise TouchableChange

  fParticleChange.ProposeTrackStatus(track.GetTrackStatus()) ;

  // If the Step was determined by the volume boundary,
  // logically relocate the particle

  if(fGeometryLimitedStep)
  {  
    // fCurrentTouchable will now become the previous touchable, 
    // and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)

    fLinearNavigator->SetGeometricallyLimitedStep() ;
    fLinearNavigator->
    LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
                                               track.GetMomentumDirection(),
                                               fCurrentTouchableHandle,
                                               true                      ) ;
    // Check whether the particle is out of the world volume 
    // If so it has exited and must be killed.
    //
    if( fCurrentTouchableHandle->GetVolume() == 0 )
    {
       fParticleChange.ProposeTrackStatus( fStopAndKill ) ;
    }
    retCurrentTouchable = fCurrentTouchableHandle ;
    fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;

    // Update the Step flag which identifies the Last Step in a volume
    if( !fFieldExertedForce )
       isLastStep =  fLinearNavigator->ExitedMotherVolume() 
                   | fLinearNavigator->EnteredDaughterVolume() ;
    else
       isLastStep = fFieldPropagator->IsLastStepInVolume(); 
  }
  else                 // fGeometryLimitedStep  is false
  {                    
    // This serves only to move the Navigator's location
    //
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

    // The value of the track's current Touchable is retained. 
    // (and it must be correct because we must use it below to
    // overwrite the (unset) one in particle change)
    //  It must be fCurrentTouchable too ??
    //
    fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
    retCurrentTouchable = track.GetTouchableHandle() ;

    isLastStep= false;
  }         // endif ( fGeometryLimitedStep ) 
  fLastStepInVolume= isLastStep; 
  
  fParticleChange.ProposeFirstStepInVolume(fFirstStepInVolume);
  fParticleChange.ProposeLastStepInVolume(isLastStep);    

  const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume() ;
  const G4Material* pNewMaterial   = 0 ;
  const G4VSensitiveDetector* pNewSensitiveDetector   = 0 ;
                                                                                       
  if( pNewVol != 0 )
  {
    pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial();
    pNewSensitiveDetector= pNewVol->GetLogicalVolume()->GetSensitiveDetector();
  }

  fParticleChange.SetMaterialInTouchable( (G4Material *) pNewMaterial ) ;
  fParticleChange.SetSensitiveDetectorInTouchable( (G4VSensitiveDetector *) pNewSensitiveDetector ) ;

  const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
  if( pNewVol != 0 )
  {
    pNewMaterialCutsCouple=pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();
  }

  if ( pNewVol!=0 && pNewMaterialCutsCouple!=0
    && pNewMaterialCutsCouple->GetMaterial()!=pNewMaterial )
  {
    // for parametrized volume
    //
    pNewMaterialCutsCouple =
      G4ProductionCutsTable::GetProductionCutsTable()
                             ->GetMaterialCutsCouple(pNewMaterial,
                               pNewMaterialCutsCouple->GetProductionCuts());
  }
  fParticleChange.SetMaterialCutsCoupleInTouchable( pNewMaterialCutsCouple );

  // temporarily until Get/Set Material of ParticleChange, 
  // and StepPoint can be made const. 
  // Set the touchable in ParticleChange
  // this must always be done because the particle change always
  // uses this value to overwrite the current touchable pointer.
  //
  fParticleChange.SetTouchableHandle(retCurrentTouchable) ;

  return &fParticleChange ;
}

/////////////////////////////////////////////////////////////////////////////
// New method takes over the responsibility to reset the state of
// G4Transportation object at the start of a new track or the resumption
// of a suspended track. 
//

void 
G4Transportation::StartTracking(G4Track* aTrack)
{
  G4VProcess::StartTracking(aTrack);
  fNewTrack= true;
  fFirstStepInVolume= true;
  fLastStepInVolume= false;
  
  // The actions here are those that were taken in AlongStepGPIL
  // when track.GetCurrentStepNumber()==1

  // reset safety value and center
  //
  fPreviousSafety    = 0.0 ; 
  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
  
  // reset looping counter -- for motion in field
  fNoLooperTrials= 0; 
  // Must clear this state .. else it depends on last track's value
  //  --> a better solution would set this from state of suspended track TODO ? 
  // Was if( aTrack->GetCurrentStepNumber()==1 ) { .. }

  // ChordFinder reset internal state
  //
  if( DoesGlobalFieldExist() )
  {
     fFieldPropagator->ClearPropagatorState();   
       // Resets all state of field propagator class (ONLY) including safety
       // values (in case of overlaps and to wipe for first track).
  }

  // Make sure to clear the chord finders of all fields (i.e. managers)
  //
  G4FieldManagerStore* fieldMgrStore = G4FieldManagerStore::GetInstance();
  fieldMgrStore->ClearAllChordFindersState(); 

  // Update the current touchable handle  (from the track's)
  //
  fCurrentTouchableHandle = aTrack->GetTouchableHandle();

  // Inform field propagator of new track
  //
  fFieldPropagator->PrepareNewTrack();
}

/////////////////////////////////////////////////////////////////////////////
//

G4bool G4Transportation::EnableUseMagneticMoment(G4bool useMoment)
{
  G4bool lastValue= fUseMagneticMoment;
  fUseMagneticMoment= useMoment;
  G4CoupledTransportation::fUseMagneticMoment= useMoment;
  return lastValue;
}

/////////////////////////////////////////////////////////////////////////////
//
void G4Transportation::SetHighLooperThresholds()
{
  // Restores the old high values -- potentially appropriate for energy-frontier
  //   HEP experiments.
  // Caution:  All tracks with E < 100 MeV that are found to loop are 
  SetThresholdWarningEnergy(    100.0 * CLHEP::MeV ); // Warn above this energy
  SetThresholdImportantEnergy(  250.0 * CLHEP::MeV ); // Extra trial above this En

  G4int maxTrials = 10;
  SetThresholdTrials( maxTrials );
  
  PushThresholdsToLogger();  // Again, to be sure
  if( verboseLevel )  ReportLooperThresholds();    
}

/////////////////////////////////////////////////////////////////////////////
void G4Transportation::SetLowLooperThresholds() // Values for low-E applications
{
  // These values were the default in Geant4 10.5 - beta
  SetThresholdWarningEnergy(     1.0 * CLHEP::keV ); // Warn above this En
  SetThresholdImportantEnergy(   1.0 * CLHEP::MeV ); // Extra trials above it

  G4int maxTrials = 30; //  A new value - was 10
  SetThresholdTrials( maxTrials );

  PushThresholdsToLogger();  // Again, to be sure
  if( verboseLevel )  ReportLooperThresholds();  
}

/////////////////////////////////////////////////////////////////////////////
//
void
G4Transportation::ReportMissingLogger( const char* methodName )
{
   const char* message= "Logger object missing from G4Transportation object";
   G4String classAndMethod= G4String("G4Transportation") + G4String( methodName );
   G4Exception(classAndMethod, "Missing Logger", JustWarning, message);   
}


/////////////////////////////////////////////////////////////////////////////
//
void
G4Transportation::ReportLooperThresholds()
{
   PushThresholdsToLogger();  // To be absolutely certain they are in sync   
   fpLogger->ReportLooperThresholds("G4Transportation");
}

/////////////////////////////////////////////////////////////////////////////
//
void G4Transportation::ProcessDescription(std::ostream& outStr) const 
 
// StreamInfo(std::ostream& out, const G4ParticleDefinition& part, G4bool rst) const
                                  
{
  G4String indent = "  "; //  : "");
  G4int oldPrec= outStr.precision(6);
  // outStr << std::setprecision(6);
  outStr << G4endl << indent << GetProcessName() << ": ";

  outStr << "   Parameters for looping particles: " << G4endl
         << "     warning-E = " << fThreshold_Warning_Energy / CLHEP::MeV  << " MeV "  << G4endl
         << "     important E = " << fThreshold_Important_Energy / CLHEP::MeV << " MeV " << G4endl
         << "     thresholdTrials " << fThresholdTrials << G4endl;
  outStr.precision(oldPrec);
}
