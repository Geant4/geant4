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
// $Id: G4CoupledTransportation.cc,v 1.28 2011-01-05 00:59:03 asaim Exp $
// --> Merged with 1.60.4.2.2.3 2007/05/09 09:30:28 japost 
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// ------------------------------------------------------------
//  GEANT 4 class implementation
// =======================================================================
// Modified:
//            13 May  2006, J. Apostolakis: Revised for parallel navigation (PathFinder)
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
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "G4CoupledTransportation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManagerStore.hh"
class G4VSensitiveDetector;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4CoupledTransportation::G4CoupledTransportation( G4int verboseLevel )
  : G4VProcess( G4String("CoupledTransportation"), fTransportation ),
    fParticleIsLooping( false ),
    fPreviousSftOrigin (0.,0.,0.),
    fPreviousMassSafety    ( 0.0 ),
    fPreviousFullSafety    ( 0.0 ),
    fThreshold_Warning_Energy( 100 * MeV ),  
    fThreshold_Important_Energy( 250 * MeV ), 
    fThresholdTrials( 10 ), 
    fUnimportant_Energy( 1 * MeV ), 
    fNoLooperTrials(0),
    fSumEnergyKilled( 0.0 ), fMaxEnergyKilled( 0.0 ), 
    fVerboseLevel(verboseLevel)     //  ( verboseLevel ? verboseLevel : 2 )  // or 4
{
  G4TransportationManager* transportMgr ; 

  transportMgr = G4TransportationManager::GetTransportationManager() ; 

  fMassNavigator = transportMgr->GetNavigatorForTracking() ; 
  fFieldPropagator = transportMgr->GetPropagatorInField() ;
  // fGlobalFieldMgr = transportMgr->GetFieldManager() ;
  fNavigatorId= transportMgr->ActivateNavigator( fMassNavigator ); 
  if( fVerboseLevel > 0 ){
    G4cout << " G4CoupledTransportation constructor: ----- " << G4endl;
    G4cout << " Verbose level is " << fVerboseLevel << G4endl;
    G4cout << " Navigator Id obtained in G4CoupledTransportation constructor " 
	   << fNavigatorId << G4endl;
  }
  fPathFinder=  G4PathFinder::GetInstance(); 
  fpSafetyHelper = transportMgr->GetSafetyHelper();  // New 

  // Following assignment is to fix small memory leak from simple use of 'new'
  static G4TouchableHandle nullTouchableHandle;  // Points to (G4VTouchable*) 0
  fCurrentTouchableHandle = nullTouchableHandle; 
  // fCurrentTouchableHandle = G4TouchableHandle( 0 );  // new G4TouchableHistory();

  fEndGlobalTimeComputed  = false;
  fCandidateEndGlobalTime = 0;
}

//////////////////////////////////////////////////////////////////////////

G4CoupledTransportation::~G4CoupledTransportation()
{
  // fCurrentTouchableHandle is a data member - no deletion required

  if( (fVerboseLevel > 0) || (fSumEnergyKilled > 0.0 ) ){ 
    G4cout << " G4CoupledTransportation: Statistics for looping particles " << G4endl;
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

G4double G4CoupledTransportation::
AlongStepGetPhysicalInteractionLength( const G4Track&  track,
                                             G4double, //  previousStepSize
                                             G4double  currentMinimumStep,
                                             G4double& proposedSafetyForStart,
                                             G4GPILSelection* selection )
{
  G4double geometryStepLength; 
  G4double startMassSafety= 0.0;   //  estimated safety for start point (mass geometry)
  G4double startFullSafety= 0.0;   //  estimated safety for start point (all geometries)
  G4double safetyProposal= -1.0;   //  local copy of proposal 

  G4ThreeVector  EndUnitMomentum ;
  G4double       lengthAlongCurve=0.0 ;
 
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
  G4VPhysicalVolume* currentVolume= track.GetVolume(); 

#ifdef G4DEBUG_TRANSPORT
  if( fVerboseLevel > 1 ) {
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
  startMassSafety = 0.0; 
  startFullSafety= 0.0; 

  //  Recall that FullSafety <= MassSafety 
  // Original: if( MagSqShift < sqr(fPreviousMassSafety) ) {
  if( MagSqShift < sqr(fPreviousFullSafety) ) {  // Revision proposed by Alex H, 2 Oct 07
     G4double mag_shift= std::sqrt(MagSqShift); 
     startMassSafety = std::max( (fPreviousMassSafety - mag_shift), 0.0); 
     startFullSafety = std::max( (fPreviousFullSafety - mag_shift), 0.0);
       // Need to be consistent between full safety with Mass safety
       //   in order reproduce results in simple case  --> use same calculation method

     // Only compute full safety if massSafety > 0.  Else it remains 0
     //   startFullSafety = fPathFinder->ComputeSafety( startPosition ); 
  }

  // Is the particle charged ?
  //
  G4double              particleCharge = pParticle->GetCharge() ; 

  fMassGeometryLimitedStep = false ; //  Set default - alex
  fAnyGeometryLimitedStep = false; 

  // fEndGlobalTimeComputed = false ;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"
  // Check whether the particle has an (EM) field force exerted upon it
  //
  G4FieldManager* fieldMgr=0;
  G4bool fieldExertsForce = false; 
  if( (particleCharge != 0.0) ) // ||  (magneticMoment != 0.0 ) )
  {
     fieldMgr= fFieldPropagator->FindAndSetFieldManager( currentVolume ); 
     if (fieldMgr != 0) {
	// Message the field Manager, to configure it for this track
	fieldMgr->ConfigureForTrack( &track );
        fieldExertsForce = (fieldMgr->GetDetectorField() != 0); 
     } 
     // the PathFinder will recognise whether the field exerts force
  }
  G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
  G4double       restMass = pParticleDef->GetPDGMass() ;
 
  fFieldPropagator->SetChargeMomentumMass( particleCharge,    // in e+ units
					   momentumMagnitude, // in Mev/c 
					   restMass           ) ;  
  // This should be obsolete - the call to SetChargeAndMoments below should do the work

  G4ThreeVector spin        = track.GetPolarization() ;
  G4FieldTrack  theFieldTrack = G4FieldTrack( startPosition, 
					    track.GetMomentumDirection(),
					    0.0, 
					    track.GetKineticEnergy(),
					    restMass,
					    0.0,    // UNUSED: track.GetVelocity(),
					    track.GetGlobalTime(), // Lab.
					    track.GetProperTime(), // Part.
					    &spin                  ) ;
  theFieldTrack.SetChargeAndMoments( particleCharge ); // EM moments -- future extension 

  G4int stepNo= track.GetCurrentStepNumber(); 

  ELimited limitedStep; 
  G4FieldTrack endTrackState('a');  //  Default values

  fMassGeometryLimitedStep = false ;    //  default 
  fAnyGeometryLimitedStep  = false ;
  if( currentMinimumStep > 0 )  {
      G4double newMassSafety= 0.0;     //  temp. for recalculation

      // Do the Transport in the field (non recti-linear)
      //
      lengthAlongCurve = fPathFinder->ComputeStep( theFieldTrack,
						   currentMinimumStep, 
						   fNavigatorId,
						   stepNo,
						   newMassSafety,
						   limitedStep,
						   endTrackState,
						   currentVolume ) ;
      // G4cout << " PathFinder ComputeStep returns " << lengthAlongCurve << G4endl; 

      G4double newFullSafety= fPathFinder->GetCurrentSafety();  
               // this was estimated already in step above
      // G4double newFullStep= fPathFinder->GetMinimumStep(); 

      if( limitedStep == kUnique || limitedStep == kSharedTransport ) {
	 fMassGeometryLimitedStep = true ;
      }

      fAnyGeometryLimitedStep = (fPathFinder->GetNumberGeometriesLimitingStep() != 0); 

      // #ifdef G4DEBUG
      if( fMassGeometryLimitedStep && !fAnyGeometryLimitedStep ){
         G4cerr << " Error in determining geometries limiting the step" << G4endl;
	 G4cerr << "  Limiting:  mass=" << fMassGeometryLimitedStep
		<< " any= " << fAnyGeometryLimitedStep << G4endl;
	 G4Exception("G4CoupledTransportation::AlongStepGetPhysicalInteractionLength()", 
		     "PathFinderConfused", 
		     FatalException, 
		     "Incompatible conditions - was limited by a geometry?");
      }
      // #endif

      // Other potential 
      // fAnyGeometryLimitedStep = newFullStep < currentMinimumStep; 
      //                                      ^^^ Not good enough; 
      //          Must compare with maximum requested step size
      //           (eg in case another process requested bigger, got this!)

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
      if( fVerboseLevel > 1 ){
	G4cout << "G4Transport:CompStep> " 
	       << " called the pathfinder for a new step at " << startPosition
	       << " and obtained step = " << lengthAlongCurve << G4endl;
	G4cout << "  New safety (preStep) = " << newMassSafety 
	       << " versus precalculated = "  << startMassSafety << G4endl; 
      }
#endif

      // Store as best estimate value
      startMassSafety    = newMassSafety ; 
      startFullSafety    = newFullSafety ; 

      // Get the End-Position and End-Momentum (Dir-ection)
      fTransportEndPosition = endTrackState.GetPosition() ;
      fTransportEndKineticEnergy  = endTrackState.GetKineticEnergy() ; 
  } else {
      geometryStepLength   = lengthAlongCurve= 0.0 ;
      fMomentumChanged         = false ; 
      // fMassGeometryLimitedStep = false ;   //  --- ???
      // fAnyGeometryLimitedStep = true;
      fTransportEndMomentumDir = track.GetMomentumDirection();
      fTransportEndKineticEnergy  = track.GetKineticEnergy();

      fTransportEndPosition = startPosition;
      // If the step length requested is 0, and we are on a boundary
      //   then a boundary will also limit the step.
      if( startMassSafety == 0.0 )  {
	 fMassGeometryLimitedStep = true ;
	 fAnyGeometryLimitedStep = true;
      }
      //   TODO:  Add explicit logical status for being at a boundary
  }
  // G4FieldTrack aTrackState(endTrackState);  

  if( !fieldExertsForce ) 
  { 
      fParticleIsLooping         = false ; 
      fMomentumChanged           = false ; 
      fEndGlobalTimeComputed     = false ; 
      // G4cout << " global time is false " << G4endl; 
  } 
  else 
  { 
  
#ifdef G4DEBUG_TRANSPORT
      if( fVerboseLevel > 1 ){
	G4cout << " G4CT::CS End Position = "  << fTransportEndPosition << G4endl; 
	G4cout << " G4CT::CS End Direction = " << fTransportEndMomentumDir << G4endl; 
      }
#endif
      //      G4cout << " G4CoupledTransportation Before if change energy statement: " << fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() << G4endl;
      if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() )
      {
	  // If the field can change energy, then the time must be integrated
	  //    - so this should have been updated
	  //
	  fCandidateEndGlobalTime   = endTrackState.GetLabTimeOfFlight();
	  fEndGlobalTimeComputed    = true;
	  //	  G4cout << " setting global time to true - why? " << G4endl;
	  
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
      
	  static G4int no_inexact_steps=0; // , no_large_ediff;
	  G4double absEdiff = std::fabs(startEnergy- endEnergy);
	  if( absEdiff > perMillion * endEnergy )  {
	    no_inexact_steps++;
	    // Possible statistics keeping here ...
	  }
#ifdef G4VERBOSE
	  if( (fVerboseLevel > 1) && ( absEdiff > perThousand * endEnergy) ){
	    ReportInexactEnergy(startEnergy, endEnergy); 
	  }  // end of if (fVerboseLevel)
#endif
	  
	  // Correct the energy for fields that conserve it
	  //  This - hides the integration error
	  //       - but gives a better physical answer
	  fTransportEndKineticEnergy= track.GetKineticEnergy(); 
      }
  }

  endpointDistance   = (fTransportEndPosition - startPosition).mag() ;
  fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;

  fTransportEndSpin = endTrackState.GetSpin();

  // Calculate the safety 
  safetyProposal= startFullSafety;   // used to be startMassSafety
     // Changed to accomodate processes that cannot update the safety -- JA 22 Nov 06

  // Update safety for the end-point, if becomes negative at the end-point.

  if(   (startFullSafety < endpointDistance ) 
	&& ( particleCharge != 0.0 ) )        //  Only needed to prepare for Mult Scat.
   //   && !fAnyGeometryLimitedStep )          // To-Try:  No safety update if at a boundary
  {
      G4double endFullSafety =
	fPathFinder->ComputeSafety( fTransportEndPosition); 
        // Expected mission -- only mass geometry's safety
        //   fMassNavigator->ComputeSafety( fTransportEndPosition) ;
        // Yet discrete processes only have poststep -- and this cannot 
        //   currently revise the safety  
        //   ==> so we use the all-geometry safety as a precaution

      fpSafetyHelper->SetCurrentSafety( endFullSafety, fTransportEndPosition);
        // Pushing safety to Helper avoids recalculation at this point

      G4ThreeVector centerPt= G4ThreeVector(0.0, 0.0, 0.0);  // Used for return value
      G4double endMassSafety= fPathFinder->ObtainSafety( fNavigatorId, centerPt); 
        //  Retrieves the mass value from PathFinder (it calculated it)

      fPreviousMassSafety = endMassSafety ; 
      fPreviousFullSafety = endFullSafety; 
      fPreviousSftOrigin = fTransportEndPosition ;

      // The convention (Stepping Manager's) is safety from the start point
      //
      safetyProposal = endFullSafety + endpointDistance;
          //  --> was endMassSafety
      // Changed to accomodate processes that cannot update the safety -- JA 22 Nov 06

      // #define G4DEBUG_TRANSPORT 1

#ifdef G4DEBUG_TRANSPORT 
      int prec= G4cout.precision(12) ;
      G4cout << "***Transportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Revised Safety at endpoint "  << fTransportEndPosition
             << "   give safety values: Mass= " << endMassSafety 
	     << "  All= " << endFullSafety << G4endl ; 
      G4cout << "  Adding endpoint distance " << endpointDistance 
             << "   to obtain pseudo-safety= " << safetyProposal << G4endl ; 
      G4cout.precision(prec); 
  }  
  else{
      int prec= G4cout.precision(12) ;
      G4cout << "***Transportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Quick Safety estimate at endpoint "  << fTransportEndPosition
	     << "   gives safety endpoint value = " << startFullSafety - endpointDistance
	     << "  using start-point value " << startFullSafety 
	     << "  and endpointDistance " << endpointDistance << G4endl; 
      G4cout.precision(prec); 
#endif
  }          

  proposedSafetyForStart= safetyProposal; 
  fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

  return geometryStepLength ;
}

//////////////////////////////////////////////////////////////////////////

G4VParticleChange* G4CoupledTransportation::AlongStepDoIt( const G4Track& track,
                                                    const G4Step&  stepData )
{
  static G4int noCalls=0;
  static const G4ParticleDefinition* fOpticalPhoton =
           G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

  noCalls++;

  fParticleChange.Initialize(track) ;
      // sets all its members to the value of corresponding members in G4Track

  //  Code specific for Transport
  //
  fParticleChange.ProposePosition(fTransportEndPosition) ;
  // G4cout << " G4CoupledTransportation::AlongStepDoIt" 
  //     << " proposes position = " << fTransportEndPosition  
  //     << " and end momentum direction  = " << fTransportEndMomentumDir <<  G4endl;
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
    G4double stepLength      = track.GetStepLength();
    
    deltaTime= 0.0;  // in case initialVelocity = 0 
    
    const G4DynamicParticle* fpDynamicParticle = track.GetDynamicParticle();
    if (fpDynamicParticle->GetDefinition()== fOpticalPhoton){
      //  A photon is in the medium of the final point
      //  during the step, so it has the final velocity.
      G4double finalVelocity   = track.GetVelocity();
      deltaTime = stepLength/finalVelocity ;
    } else {
      G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity();
      if ( initialVelocity > 0.0 ) deltaTime = stepLength/initialVelocity ;
    }
    fCandidateEndGlobalTime   = startTime + deltaTime ;
  }
  else
  {
     deltaTime = fCandidateEndGlobalTime - startTime ;
     // G4cout << " Calculated global time from candidate end time = "
     //    << fCandidateEndGlobalTime << " and start time = " << startTime << G4endl;
  }

  // G4cout << " G4CoupledTransportation::AlongStepDoIt  "
  // << " flag whether computed time = " << fEndGlobalTimeComputed  << " and " 
  // << " is proposes end time " << fCandidateEndGlobalTime << G4endl; 
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
	  G4cout << " G4CoupledTransportation is killing track that is looping or stuck " << G4endl
		 << "   This track has " << track.GetKineticEnergy() / MeV
		 << " MeV energy." << G4endl;
	}
	if( fVerboseLevel > 0 ) { 
	  G4cout << "   Steps by this track: " << track.GetCurrentStepNumber() << G4endl;
	}
#endif
	fNoLooperTrials=0; 
      } else{ 
	fNoLooperTrials ++; 
#ifdef G4VERBOSE
	if( (fVerboseLevel > 2) ){
	  G4cout << "  ** G4CoupledTransportation::AlongStepDoIt(): Particle looping -  " << G4endl
		 << "   Number of consecutive problem step (this track) = " << fNoLooperTrials << G4endl
		 << "   Steps by this track: " << track.GetCurrentStepNumber() << G4endl
		 << "   Total no of calls to this method (all tracks) = " << noCalls << G4endl;
	}
#endif
      }
  }else{ 
      fNoLooperTrials=0; 
  }

  // Another (sometimes better way) is to use a user-limit maximum Step size
  // to alleviate this problem .. 

  // Add smooth curved trajectories to particle-change
  //
  // fParticleChange.SetPointerToVectorOfAuxiliaryPoints
  //   (fFieldPropagator->GimmeTrajectoryVectorAndForgetIt() );

  return &fParticleChange ;
}

//////////////////////////////////////////////////////////////////////////
//
//  This ensures that the PostStep action is always called,
//  so that it can do the relocation if it is needed.
// 

G4double G4CoupledTransportation::
PostStepGetPhysicalInteractionLength( const G4Track&,
                                            G4double, // previousStepSize
                                            G4ForceCondition* pForceCond )
{ 
  // Must act as PostStep action -- to relocate particle
  *pForceCond = Forced ;    
  return DBL_MAX ;
}

static 
void ReportMove( G4ThreeVector OldVector, G4ThreeVector NewVector, G4String& Quantity )
{
    G4ThreeVector moveVec = ( NewVector - OldVector );

    G4cerr << G4endl
	   << "**************************************************************" << G4endl;
    G4cerr << "Endpoint has moved between value expected from TransportEndPosition "
	   << " and value from Track in PostStepDoIt. " << G4endl
	   << "Change of " << Quantity << " is " << moveVec.mag() / mm << " mm long, "
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

  // Check that the end position and direction are preserved 
  //   since call to AlongStepDoIt
  if( (fTransportEndPosition  - track.GetPosition()).mag2() >= 1.0e-16 ){
     static G4String EndLabelString("End of Step Position");  
     ReportMove( track.GetPosition(), fTransportEndPosition, EndLabelString ); 
     G4cerr << " Problem in G4CoupledTransportation::PostStepDoIt " << G4endl; 
  }

  // If the Step was determined by the volume boundary, relocate the particle
  // The pathFinder will know that the geometry limited the step (!?)

  if( fVerboseLevel > 0 ){
     G4cout << " Calling PathFinder::Locate() from " 
	    << " G4CoupledTransportation::PostStepDoIt " << G4endl;
     G4cout << "  fAnyGeometryLimitedStep is " << fAnyGeometryLimitedStep << G4endl;

  }
  if(fAnyGeometryLimitedStep)
  {  
    fPathFinder->Locate( track.GetPosition(), 
		       track.GetMomentumDirection(),
		       true); 

    // fCurrentTouchable will now become the previous touchable, 
    // and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)
    if( fVerboseLevel > 0 )
      G4cout << "G4CoupledTransportation::PostStepDoIt --- fNavigatorId = " 
	     << fNavigatorId << G4endl;

    fCurrentTouchableHandle= 
      fPathFinder->CreateTouchableHandle( fNavigatorId );

    // Check whether the particle is out of the world volume 
    // If so it has exited and must be killed.
    //
#ifdef G4DEBUG_TRANSPORT
    if( fVerboseLevel > 1 ){
       G4VPhysicalVolume* vol= fCurrentTouchableHandle->GetVolume(); 
       G4cout << "CHECK !!!!!!!!!!! fCurrentTouchableHandle->GetVolume() = " << vol;
       if( vol ) { G4cout << "Name=" << vol->GetName(); }
       G4cout << G4endl;
    }
#endif
    if( fCurrentTouchableHandle->GetVolume() == 0 )
    {
       fParticleChange.ProposeTrackStatus( fStopAndKill ) ;
    }
    retCurrentTouchable = fCurrentTouchableHandle ;
    // fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;

    // Notify particle change that this is last step in volume
    fParticleChange.ProposeLastStepInVolume(true);
    // Double check that a boundary limited the step, and 
    // if( fLinearNavigator->Get

  }
  else                 // fAnyGeometryLimitedStep  is false
  { 
#ifdef G4DEBUG_TRANSPORT
    if( fVerboseLevel > 1 ){
       G4cout << "G4CoupledTransportation::PostStepDoIt -- "
              << " fAnyGeometryLimitedStep  = " << fAnyGeometryLimitedStep  
              << " must be false " << G4endl;
    }
#endif
    // This serves only to move each of the Navigator's location
    //
    // fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

    // G4cout << "G4CoupledTransportation calling PathFinder::ReLocate() " << G4endl;
    fPathFinder->ReLocate( track.GetPosition() );
			   // track.GetMomentumDirection() ); 

    // Keep the value of the track's current Touchable is retained,
    //  and use it to overwrite the (unset) one in particle change.
    // Expect this must be fCurrentTouchable too
    //   - could it be different, eg at the start of a step ?
    //
    retCurrentTouchable = track.GetTouchableHandle() ;
    // fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;

    // Have not reached a boundary
    fParticleChange.ProposeLastStepInVolume(false);
  }         // endif ( fAnyGeometryLimitedStep ) 

  const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume() ;
  const G4Material* pNewMaterial   = 0 ;
  const G4VSensitiveDetector* pNewSensitiveDetector   = 0 ;
                                                                                       
  if( pNewVol != 0 )
  {
    pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial();
    pNewSensitiveDetector= pNewVol->GetLogicalVolume()->GetSensitiveDetector();
  }

  // ( const_cast<G4Material *> pNewMaterial ) ;
  // ( const_cast<G4VSensitiveDetetor *> pNewSensitiveDetector) ;

  fParticleChange.SetMaterialInTouchable( (G4Material *) pNewMaterial ) ;
  fParticleChange.SetSensitiveDetectorInTouchable( (G4VSensitiveDetector *) pNewSensitiveDetector ) ;
             // "temporarily" until Get/Set Material of ParticleChange, 
             // and StepPoint can be made const. 

  const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
  if( pNewVol != 0 )
  {
    pNewMaterialCutsCouple=pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();
    if( pNewMaterialCutsCouple!=0 
	&& pNewMaterialCutsCouple->GetMaterial()!=pNewMaterial )
      {
	// for parametrized volume
	//
	pNewMaterialCutsCouple =
	  G4ProductionCutsTable::GetProductionCutsTable()
                       ->GetMaterialCutsCouple(pNewMaterial,
					       pNewMaterialCutsCouple->GetProductionCuts());
      }
  }
  fParticleChange.SetMaterialCutsCoupleInTouchable( pNewMaterialCutsCouple );

  // Must always set the touchable in ParticleChange, whether relocated or not
  fParticleChange.SetTouchableHandle(retCurrentTouchable) ;

  return &fParticleChange ;
}

// New method takes over the responsibility to reset the state of 
//   G4CoupledTransportation object:
//      - at the start of a new track,  and
//      - on the resumption of a suspended track. 

void 
G4CoupledTransportation::StartTracking(G4Track* aTrack)
{

  static G4TransportationManager* transportMgr=  
      G4TransportationManager::GetTransportationManager(); 

  // G4VProcess::StartTracking(aTrack);

  //  The 'initialising' actions
  //     once taken in AlongStepGPIL -- if ( track.GetCurrentStepNumber()==1 )

  // fStartedNewTrack= true; 

  fMassNavigator = transportMgr->GetNavigatorForTracking() ; 
  fNavigatorId= transportMgr->ActivateNavigator( fMassNavigator );  // Confirm it!

  // if( fVerboseLevel > 1 ){
  //  G4cout << " Navigator Id obtained in StartTracking " << fNavigatorId << G4endl;
  // }
  G4ThreeVector position = aTrack->GetPosition(); 
  G4ThreeVector direction = aTrack->GetMomentumDirection();

  // if( fVerboseLevel > 1 ){
  //   G4cout << " Calling PathFinder::PrepareNewTrack from    " 
  //   << " G4CoupledTransportation::StartTracking -- which calls Locate()" << G4endl;
  // }
  fPathFinder->PrepareNewTrack( position, direction); 
  // This implies a call to fPathFinder->Locate( position, direction ); 

  // Global field, if any, must exist before tracking is started
  fGlobalFieldExists= DoesGlobalFieldExist(); 
  // reset safety value and center
  //
  fPreviousMassSafety  = 0.0 ; 
  fPreviousFullSafety  = 0.0 ; 
  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
  
  // reset looping counter -- for motion in field  
  fNoLooperTrials= 0; 
  // Must clear this state .. else it depends on last track's value
  //  --> a better solution would set this from state of suspended track TODO ? 
  // Was if( aTrack->GetCurrentStepNumber()==1 ) { .. }

  // ChordFinder reset internal state
  //
  if( fGlobalFieldExists ) {
     fFieldPropagator->ClearPropagatorState();   
       // Resets safety values, in case of overlaps.  

     G4ChordFinder* chordF= fFieldPropagator->GetChordFinder();
     if( chordF ) chordF->ResetStepEstimate();
  }
  // Clear the chord finders of all fields (ie managers) derived objects
  static G4FieldManagerStore* fieldMgrStore= G4FieldManagerStore::GetInstance();
  fieldMgrStore->ClearAllChordFindersState(); 

#ifdef G4DEBUG_TRANSPORT
  if( fVerboseLevel > 1 ){
    G4cout << " Returning touchable handle " << fCurrentTouchableHandle << G4endl;
  }
#endif

  // Update the current touchable handle  (from the track's)
  //
  fCurrentTouchableHandle = aTrack->GetTouchableHandle();  
}

void 
G4CoupledTransportation::EndTracking()
{
  G4TransportationManager::GetTransportationManager()->InactivateAll();
}

void
G4CoupledTransportation::
ReportInexactEnergy(G4double startEnergy, G4double endEnergy)
{
  static G4int no_warnings= 0, warnModulo=1,  moduloFactor= 10, no_large_ediff= 0; 

  if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
    {
      no_large_ediff ++;
      if( (no_large_ediff% warnModulo) == 0 )
	{
	  no_warnings++;
	  G4cout << "WARNING - G4CoupledTransportation::AlongStepGetPIL() " 
		 << "   Energy change in Step is above 1^-3 relative value. " << G4endl
		 << "   Relative change in 'tracking' step = " 
		 << std::setw(15) << (endEnergy-startEnergy)/startEnergy << G4endl
		 << "     Starting E= " << std::setw(12) << startEnergy / MeV << " MeV " << G4endl
		 << "     Ending   E= " << std::setw(12) << endEnergy   / MeV << " MeV " << G4endl;       
	  G4cout << " Energy has been corrected -- however, review"
		 << " field propagation parameters for accuracy."  << G4endl;
	  if( (fVerboseLevel > 2 ) || (no_warnings<4) || (no_large_ediff == warnModulo * moduloFactor) ){
	    G4cout << " These include EpsilonStepMax(/Min) in G4FieldManager "
		   << " which determine fractional error per step for integrated quantities. " << G4endl
		   << " Note also the influence of the permitted number of integration steps."
		   << G4endl;
	  }
	  G4cerr << "ERROR - G4CoupledTransportation::AlongStepGetPIL()" << G4endl
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
}
