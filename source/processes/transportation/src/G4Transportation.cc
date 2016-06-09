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
// $Id: G4Transportation.cc,v 1.55.2.1 2006/06/29 21:12:28 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
// It is also tasked with part of updating the "safety".
//
// =======================================================================
// Modified:   
//            19 Jan  2006, P.MoraDeFreitas: Fix for suspended tracks (StartTracking)
//            11 Aug  2004, M.Asai: Add G4VSensitiveDetector* for updating stepPoint.
//            21 June 2003, J.Apostolakis: Calling field manager with 
//                            track, to enable it to configure its accuracy
//            13 May  2003, J.Apostolakis: Zero field areas now taken into
//                            account correclty in all cases (thanks to W Pokorski).
//            29 June 2001, J.Apostolakis, D.Cote-Ahern, P.Gumplinger: 
//                          correction for spin tracking   
//            20 Febr 2001, J.Apostolakis:  update for new FieldTrack
//            22 Sept 2000, V.Grichine:     update of Kinetic Energy
//             9 June 1999, J.Apostolakis & S.Giani: protect full relocation
//                          in DEBUG for track  starting on surface that
//                          goes step < tolerance.
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "G4Transportation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ChordFinder.hh"
class G4VSensitiveDetector;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4Transportation::G4Transportation( G4int verboseLevel )
  : G4VProcess( G4String("Transportation"), fTransportation ),
    fParticleIsLooping( false ),
    fPreviousSftOrigin (0.,0.,0.),
    fPreviousSafety    ( 0.0 ),
    fThreshold_Warning_Energy( 100 * MeV ),  
    fThreshold_Important_Energy( 250 * MeV ), 
    fThresholdTrials( 10 ), 
    fUnimportant_Energy( 1 * MeV ), 
    fNoLooperTrials(0),
    fSumEnergyKilled( 0.0 ), fMaxEnergyKilled( 0.0 ), 
    fVerboseLevel( verboseLevel )
{
  G4TransportationManager* transportMgr ; 

  transportMgr = G4TransportationManager::GetTransportationManager() ; 

  fLinearNavigator = transportMgr->GetNavigatorForTracking() ; 

  // fGlobalFieldMgr = transportMgr->GetFieldManager() ;

  fFieldPropagator = transportMgr->GetPropagatorInField() ;

  // Cannot determine whether a field exists here,
  //  because it would only work if the field manager has informed 
  //  about the detector's field before this transportation process 
  //  is constructed.
  // Instead later the method DoesGlobalFieldExist() is called

  fCurrentTouchableHandle = new G4TouchableHistory();
  
  fEndGlobalTimeComputed  = false;
  fCandidateEndGlobalTime = 0;
}

//////////////////////////////////////////////////////////////////////////

G4Transportation::~G4Transportation()
{
  if( (fVerboseLevel > 0) && (fSumEnergyKilled > 0.0 ) ){ 
    G4cout << " G4Transportation: Statistics for looping particles " << G4endl;
    G4cout << "   Sum of energy of loopers killed: " <<  fSumEnergyKilled << G4endl;
    G4cout << "   Max energy of loopers killed: " <<  fMaxEnergyKilled << G4endl;
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
  G4double geometryStepLength, newSafety ; 
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

  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
  const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
  G4ThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  G4ThreeVector startPosition          = track.GetPosition() ;

  // G4double   theTime        = track.GetGlobalTime() ;

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

  // Is the particle charged ?
  //
  G4double              particleCharge = pParticle->GetCharge() ; 

  fGeometryLimitedStep = false ;
  // fEndGlobalTimeComputed = false ;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"

  // Check whether the particle have an (EM) field force exerting upon it
  //
  G4FieldManager* fieldMgr=0;
  G4bool          fieldExertsForce = false ;
  if( (particleCharge != 0.0) )
  {
     fieldMgr= fFieldPropagator->FindAndSetFieldManager( track.GetVolume() ); 
     if (fieldMgr != 0) {
        // If the field manager has no field, there is no field !
        fieldExertsForce = (fieldMgr->GetDetectorField() != 0);
     } 
  }

  // Choose the calculation of the transportation: Field or not 
  //
  if( !fieldExertsForce ) 
  {
     G4double linearStepLength ;
     if( currentMinimumStep <= currentSafety )
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

       // The safety at the initial point has been re-calculated:
       //
       currentSafety = newSafety ;
          
       if( linearStepLength <= currentMinimumStep)
       {
         // The geometry limits the Step size (an intersection was found.)
         //
         geometryStepLength   = linearStepLength ;
         fGeometryLimitedStep = true ;
       }
       else
       {
         // The full Step is taken.
         //
         geometryStepLength   = currentMinimumStep ;
         fGeometryLimitedStep = false ;
       }
     }
     endpointDistance = geometryStepLength ;

     // Calculate final position
     //
     fTransportEndPosition = startPosition+geometryStepLength*startMomentumDir ;

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
     G4double       restMass = pParticleDef->GetPDGMass() ;
 
     fFieldPropagator->SetChargeMomentumMass( particleCharge,    // in e+ units
                                              momentumMagnitude, // in Mev/c 
                                              restMass           ) ;  

     // Message the field Manager, to configure it for this track
     fieldMgr->ConfigureForTrack( &track );

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
        if( lengthAlongCurve < currentMinimumStep)
        {
           geometryStepLength   = lengthAlongCurve ;
           fGeometryLimitedStep = true ;
        }
        else
        {
           geometryStepLength   = currentMinimumStep ;
           fGeometryLimitedStep = false ;
        }
     }
     else
     {
        geometryStepLength   = lengthAlongCurve= 0.0 ;
        fGeometryLimitedStep = false ;
     }

     // Remember last safety origin & value.
     //
     fPreviousSftOrigin = startPosition ;
     fPreviousSafety    = currentSafety ;         
        
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

        static G4int no_inexact_steps=0, no_large_ediff;
        G4double absEdiff = std::fabs(startEnergy- endEnergy);
        if( absEdiff > perMillion * endEnergy )
        {
          no_inexact_steps++;
          // Possible statistics keeping here ...
        }
        if( fVerboseLevel > 1 )
        {
          if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
          {
            static G4int no_warnings= 0, warnModulo=1,  moduloFactor= 10; 
            no_large_ediff ++;
            if( (no_large_ediff% warnModulo) == 0 )
            {
               no_warnings++;
               G4cout << "WARNING - G4Transportation::AlongStepGetPIL()" << G4endl
	              << "   Energy changed in Step, more than 1/1000: " << G4endl
                      << "          Start= " << startEnergy   << G4endl
                      << "          End= "   << endEnergy     << G4endl
                      << "          Relative change= "
                      << (startEnergy-endEnergy)/startEnergy << G4endl;
               G4cout << " Energy has been corrected -- however, review"
                      << " field propagation parameters for accuracy." << G4endl;
               G4cerr << "ERROR - G4Transportation::AlongStepGetPIL()" << G4endl
	              << "        Bad 'endpoint'. Energy change detected"
                      << " and corrected,"                      << G4endl
                      << "        occurred already "
                      << no_large_ediff << " times." << G4endl;
               if( no_large_ediff == warnModulo * moduloFactor )
               {
                  warnModulo *= moduloFactor;
               }
            }
          }
        }  // end of if (fVerboseLevel)

        // Correct the energy for fields that conserve it
        //  This - hides the integration error
        //       - but gives a better physical answer
        fTransportEndKineticEnergy= track.GetKineticEnergy(); 
     }

     fTransportEndSpin = aFieldTrack.GetSpin();
     fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;
     endpointDistance   = (fTransportEndPosition - startPosition).mag() ;
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  // then a boundary will also limit the step -> we must flag this.
  //
  if( currentMinimumStep == 0.0 ) 
  {
      if( currentSafety == 0.0 )  fGeometryLimitedStep = true ;
  }

  // Update the safety starting from the end-point,
  // if it will become negative at the end-point.
  //
  if( currentSafety < endpointDistance ) 
  {
      G4double endSafety =
               fLinearNavigator->ComputeSafety( fTransportEndPosition) ;
      currentSafety      = endSafety ;
      fPreviousSftOrigin = fTransportEndPosition ;
      fPreviousSafety    = currentSafety ; 

      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength
      //
      currentSafety     += endpointDistance ;

#ifdef G4DEBUG_TRANSPORT 
      G4cout.precision(16) ;
      G4cout << "***Transportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Called Navigator->ComputeSafety at "
             << fTransportEndPosition
             << "    and it returned safety= " << endSafety << G4endl ; 
      G4cout << "  Adding endpoint distance " << endpointDistance 
             << "   we obtain pseudo-safety= " << currentSafety << G4endl ; 
#endif
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
     G4double finalVelocity   = track.GetVelocity() ;
     G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity() ;
     G4double stepLength      = track.GetStepLength() ;

     static const G4ParticleDefinition* fOpticalPhoton =
           G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
     const G4DynamicParticle* fpDynamicParticle = track.GetDynamicParticle();
     if (fpDynamicParticle->GetDefinition()== fOpticalPhoton)
     {
        //  A photon is in the medium of the final point
        //  during the step, so it has the final velocity.
        deltaTime = stepLength/finalVelocity ;
     }
     else if (finalVelocity > 0.0)
     {
        G4double meanInverseVelocity ;
        // deltaTime = stepLength/finalVelocity ;
        meanInverseVelocity = 0.5
                            * ( 1.0 / initialVelocity + 1.0 / finalVelocity ) ;
        deltaTime = stepLength * meanInverseVelocity ;
     }
     else
     {
        deltaTime = stepLength/initialVelocity ;
     }
     fCandidateEndGlobalTime   = startTime + deltaTime ;
  }
  else
  {
     deltaTime = fCandidateEndGlobalTime - startTime ;
  }

  fParticleChange.ProposeGlobalTime( fCandidateEndGlobalTime ) ;

  // Now Correct by Lorentz factor to get "proper" deltaTime
  
  G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
  G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;

  fParticleChange.ProposeProperTime(track.GetProperTime() + deltaProperTime) ;
  //fParticleChange. ProposeTrueStepLength( track.GetStepLength() ) ;

  // If the particle is caught looping or is stuck (in very difficult
  // boundaries) in a magnetic field (doing many steps) 
  //   THEN this kills it ...
  //
  if ( fParticleIsLooping )
  {
      G4double endEnergy= fTransportEndKineticEnergy;

      if( (endEnergy < fThreshold_Important_Energy) 
	  || (fNoLooperTrials >= fThresholdTrials ) ){
	// Kill the looping particle 
	//
	fParticleChange.ProposeTrackStatus( fStopAndKill )  ;

        // 'Bare' statistics
        fSumEnergyKilled += endEnergy; 
	if( endEnergy > fMaxEnergyKilled) { fMaxEnergyKilled= endEnergy; }

#ifdef G4VERBOSE
	if( (fVerboseLevel > 1) || 
	    ( endEnergy > fThreshold_Warning_Energy )  ) { 
	  G4cout << " G4Transportation is killing track that is looping or stuck "
		 << G4endl
		 << "   This track has " << track.GetKineticEnergy() / MeV
		 << " MeV energy." << G4endl;
	}
#endif
	fNoLooperTrials=0; 
      }
      else{
	fNoLooperTrials ++; 
      }
  }else{
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
  *pForceCond = Forced ; 
  return DBL_MAX ;  // was kInfinity ; but convention now is DBL_MAX
}

/////////////////////////////////////////////////////////////////////////////
//

G4VParticleChange* G4Transportation::PostStepDoIt( const G4Track& track,
                                                   const G4Step& )
{
  G4TouchableHandle retCurrentTouchable ;   // The one to return

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
  }
  else                 // fGeometryLimitedStep  is false
  {                    
#ifdef G4DEBUG_POSTSTEP_TRANSPORT

    // Although the location is changed, we know that the physical 
    // volume remains constant. 
    // In order to help in checking the user geometry
    // we perform a full-relocation and check its result 
    // *except* if we have made a very small step from a boundary
    // (i.e. remaining inside the tolerance)

    G4bool  startAtSurface_And_MoveEpsilon ;
    startAtSurface_And_MoveEpsilon =
             (stepData.GetPreStepPoint()->GetSafety() == 0.0) && 
             (stepData.GetStepLength() < kCarTolerance ) ;

    if( startAtSurface_And_MoveEpsilon ) 
    {
       fLinearNavigator->
       LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
                                                  track.GetMomentumDirection(),
                                                  fCurrentTouchableHandle,
                                                  true                     );
       if( fCurrentTouchableHandle->GetVolume() != track.GetVolume() )
       {
         G4cerr << " ERROR: A relocation within safety has"
                << " caused a volume change! " << G4endl  ; 
         G4cerr << "   The old volume is called " 
                << track.GetVolume()->GetName() << G4endl ; 
         G4cerr << "   The new volume is called " ;

         if ( fCurrentTouchableHandle->GetVolume() != 0 )
         {
            G4cerr << fCurrentTouchableHandle->GetVolume()->GetName()
                   << G4endl ; 
         }
         else
         {
            G4cerr << "Out of World" << G4endl ; 
         }
         G4cerr.precision(7) ;
         G4cerr << "   The position is " << track.GetPosition() <<  G4endl ;

         // Let us relocate again, for debuging
         //
         fLinearNavigator->
         LocateGlobalPointAndUpdateTouchableHandle(track.GetPosition(),
                                                   track.GetMomentumDirection(),
                                                   fCurrentTouchableHandle,
                                                   true                     ) ;
         G4cerr << "   The newer volume is called "  ;

         if ( fCurrentTouchableHandle->GetVolume() != 0 )
         {
            G4cerr << fCurrentTouchableHandle->GetVolume()->GetName()
                   << G4endl ;
         } 
         else
         {
            G4cerr << "Out of World" << G4endl ; 
         }
       }

       assert( fCurrentTouchableHandle->GetVolume()->GetName() ==
               track.GetVolume()->GetName() ) ;

       retCurrentTouchable = fCurrentTouchableHandle ; 
       fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;
       
    }
    else
    {
       retCurrentTouchable = track.GetTouchableHandle() ;
       fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
    }

    //  This must be done in the above if ( AtSur ) fails
    //  We also do it for if (true) in order to get debug/opt to  
    //  behave as exactly the same way as possible.
    //
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

#else    // ie #ifndef G4DEBUG_POSTSTEP_TRANSPORT does a quick relocation

    // This serves only to move the Navigator's location
    //
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

    // The value of the track's current Touchable is retained. 
    // (and it must be correct because we must use it below to
    // overwrite the (unset) one in particle change)
    // Although in general this is fCurrentTouchable, at the start of
    // a step it could be different ... ??
    //
    fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
    retCurrentTouchable = track.GetTouchableHandle() ;

#endif
  }         // endif ( fGeometryLimitedStep ) 

  const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume() ;
  const G4Material* pNewMaterial   = 0 ;
  const G4VSensitiveDetector* pNewSensitiveDetector   = 0 ;
                                                                                       
  if( pNewVol != 0 )
  {
    pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial();
    pNewSensitiveDetector= pNewVol->GetLogicalVolume()->GetSensitiveDetector();
  }

  // ( <const_cast> pNewMaterial ) ;
  // ( <const_cast> pNewSensitiveDetector) ;

  fParticleChange.SetMaterialInTouchable( (G4Material *) pNewMaterial ) ;
  fParticleChange.SetSensitiveDetectorInTouchable( (G4VSensitiveDetector *) pNewSensitiveDetector ) ;

  const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
  if( pNewVol != 0 )
  {
    pNewMaterialCutsCouple=pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();
  }

  if( pNewVol!=0 && pNewMaterialCutsCouple!=0 && pNewMaterialCutsCouple->GetMaterial()!=pNewMaterial )
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

// New method takes over the responsibility to reset the state of G4Transportation
//   object at the start of a new track or the resumption of a suspended track. 

void 
G4Transportation::StartTracking(G4Track* aTrack)
{
  G4VProcess::StartTracking(aTrack);

// The actions here are those that were taken in AlongStepGPIL
//   when track.GetCurrentStepNumber()==1

  // reset safety value and center
  //
  fPreviousSafety    = 0.0 ; 
  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
  
  // reset looping counter -- for motion in field
  if( aTrack->GetCurrentStepNumber()==1 ) {
     fNoLooperTrials= 0; 
  }

  // ChordFinder reset internal state
  //
  if( DoesGlobalFieldExist() ) {
     G4ChordFinder* chordF= fFieldPropagator->GetChordFinder();
     if( chordF ) chordF->ResetStepEstimate();
  }
  
  // Update the current touchable handle  (from the track's)
  //
  fCurrentTouchableHandle = aTrack->GetTouchableHandle();
}

