// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Transportation.cc,v 1.11 2000-06-19 16:13:48 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//
// 
// ------------------------------------------------------------
//	GEANT 4  include file implementation
//
//	For information related to this code contact:
//	CERN, IT Division (formely CN), ASD group
// ------------------------------------------------------------
//
//   This class is a process responsible for the transportation of 
// a particle, ie the geometrical propagation that encounters the 
// geometrical sub-volumes of the detectors.
//
//   It is also tasked with part of updating the "safety".
//
// =======================================================================
// Created:   19 March 1997, J. Apostolakis
// Modified:   9 June  1999, J. Apostolakis & S.Giani: protect full relocation used in DEBUG 
//                 		 for track that started on surface and went step < tolerance
//    				Also forced fast relocation in all DEBUG cases
//    				 & changed #if to use DEBUG instead of VERBOSE
// =======================================================================

#include "G4Transportation.hh"

G4Transportation::G4Transportation() :
      G4VProcess(G4String("Transportation") )
{
  G4TransportationManager* transportMgr;

  transportMgr= G4TransportationManager::GetTransportationManager();

  fLinearNavigator=   transportMgr->GetNavigatorForTracking();
  fFieldPropagator=   0;

  // fFieldExists= false;
  fParticleIsLooping = false;
 
  // fGlobalFieldMgr=    transportMgr->GetFieldManager();
  fFieldPropagator=   transportMgr->GetPropagatorInField();

  // Find out if an electromagnetic field exists
  // 
  // fFieldExists= transportMgr->GetFieldManager()->DoesFieldExist();
  // 
  //   The above code is problematic, because it only works if
  // the field manager has informed about the detector's field 
  // before this transportation process is constructed.
  // I cannot foresee how the transportation can be informed later. JA 
  //   The current answer is to ignore this data member and use 
  // the member function DoesGlobalFieldExist() in its place ...
  //    John Apostolakis, July 7, 1997

  fTouchable1 = new G4TouchableHistory();
  fTouchable2 = new G4TouchableHistory();

  fIsTouchable1Free= true;
  fIsTouchable2Free= true;

  // Initial value for safety and point-of-origin of safety
  fPreviousSafety=0.0; 
  fPreviousSftOrigin= G4ThreeVector(0.,0.,0.);

}

G4Transportation::~G4Transportation()
{
   delete fTouchable1;
   delete fTouchable2;
}

// ------------------------------------------------------------------
// G4double G4Transportation::GetContinuousStepLimit  (
G4double G4Transportation::AlongStepGetPhysicalInteractionLength(
	   const G4Track&  track,
		 G4double  previousStepSize,
		 G4double  currentMinimumStep,
		 G4double& currentSafety,
		 G4GPILSelection* selection
	 )
// ------------------------------------------------------------------
{
  // Responsibilities:
  //    Find whether the geometry limits the Step, and to what length
  //    Calculate the new value of the safety and return it.
  //    Store the final time, position and momentum.
  G4double geometryStepLength, newSafety; 
  fParticleIsLooping = false;

  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  *selection = CandidateForSelection;

  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle*  pParticle  = track.GetDynamicParticle();
  G4double      startEnergy      = pParticle->GetKineticEnergy();
  G4ThreeVector startMomentumDir = pParticle->GetMomentumDirection();
  G4ThreeVector startPosition  = track.GetPosition();
  // G4double   theTime        = track.GetGlobalTime();

  // The Step Point safety is now generalised to mean the limit of assumption
  // of all processes, so it is not the previous Step's geometrical safety.
  //
  // We calculate the starting point's safety here.
  G4ThreeVector OriginShift= startPosition - fPreviousSftOrigin;
  G4double      MagSqShift=  OriginShift.mag2();
  if( MagSqShift >= sqr(fPreviousSafety) ){
     currentSafety = 0.0;
  }else{
     currentSafety = fPreviousSafety - sqrt(MagSqShift);
  }

  // Is the particle charged ?
  G4ParticleDefinition* pParticleDef=   pParticle->GetDefinition();
  G4double              particleCharge= pParticleDef->GetPDGCharge(); 

  G4bool   fieldExertsForce= false;
  fGeometryLimitedStep= false;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"
  //

  //  Does the particle have an (EM) field force exerting upon it?
  //
  if( (particleCharge!=0.0) ){
     
     fieldExertsForce= this->DoesGlobalFieldExist();
     // Future: will/can also check whether current volume's field is Zero or
     //  set by the user (in the logical volume) to be zero.
  }

  //  Choose the calculation of the transportation: Field or not 
  //
  if( !fieldExertsForce ) 
  {
     G4double linearStepLength;
     
     if( currentMinimumStep <= currentSafety )
     {
        // The Step is guaranteed to be taken
	geometryStepLength=currentMinimumStep;
	fGeometryLimitedStep= false;
     }
     else
     {
	//  Find whether the straight path intersects a volume
	linearStepLength= fLinearNavigator->ComputeStep( 
					 startPosition, startMomentumDir,
					 currentMinimumStep, newSafety);

        // Remember last safety origin & value.
	fPreviousSftOrigin = startPosition;
        fPreviousSafety= newSafety; 
	// The safety at the initial point has been re-calculated:
        currentSafety= newSafety;
			    
	if( linearStepLength <= currentMinimumStep){
	   // The geometry limits the Step size (an intersection was found.)
	   geometryStepLength=linearStepLength;
	   fGeometryLimitedStep= true;
	}else{
	   // The full Step is taken.
	   geometryStepLength=currentMinimumStep;
	   fGeometryLimitedStep= false;
	}
     }
     endpointDistance= geometryStepLength;

     // Calculate final position
     fTransportEndPosition= startPosition+geometryStepLength*startMomentumDir;
     // Momentum (& its direction) is unchanged
     fTransportEndMomentumDir= startMomentumDir; 
     fTransportEndKineticEnergy= track.GetKineticEnergy();
     fParticleIsLooping = false;
     fMomentumChanged = false; 
  }
  else
  {
     G4double       momentumMagnitude=pParticle->GetTotalMomentum();
     G4ThreeVector  EndUnitMomentum;
     G4double       lengthAlongCurve;
     G4double       restMass= pParticleDef->GetPDGMass();
 
     fFieldPropagator->SetChargeMomentumMass(
		     particleCharge,           // charge in e+ units
		     momentumMagnitude,        // Momentum in Mev/c 
		     restMass );  

     G4ThreeVector spin = track.GetPolarization();   // Does it have it ?
     G4ThreeVector velocityVector =   track.GetVelocity() 
                                    * track.GetMomentumDirection(); 
     G4FieldTrack  aFieldTrack = 
       G4FieldTrack(  startPosition, 
		      velocityVector,
		      0.0, 
		      track.GetKineticEnergy(),
		      track.GetLocalTime(),    // tof lab ?
		      track.GetProperTime(),   // tof proper
		      &spin );

     if( currentMinimumStep > 0 ) {
        //  Do the Transport in the field (non recti-linear)
        lengthAlongCurve=fFieldPropagator->ComputeStep( aFieldTrack,
							currentMinimumStep, 
							currentSafety,
							track.GetVolume() );
	//               ----------------
	if( lengthAlongCurve< currentMinimumStep){
	   geometryStepLength=lengthAlongCurve;
	   fGeometryLimitedStep= true;
	}else{
  	   geometryStepLength=currentMinimumStep;
	   fGeometryLimitedStep= false;
	}
     }else{
        geometryStepLength= lengthAlongCurve= 0.0;
	fGeometryLimitedStep= false;
     }

     // Remember last safety origin & value.
     fPreviousSftOrigin = startPosition;
     fPreviousSafety= currentSafety; 		    
        
     // Get the End-Position and End-Momentum (Dir-ection)
     fTransportEndPosition= aFieldTrack.GetPosition();
     // Momentum:  Magnitude and direction can be changed too now ...
     fMomentumChanged = true; 
     fTransportEndMomentumDir= aFieldTrack.GetMomentumDir();

     // fTransportEndKineticEnergy= aFieldTrack.GetEnergy(); // Energy is wrong
#if 0
     G4ThreeVector endVelocity = aFieldTrack.GetVelocity();
     G4double  veloc_sq = endVelocity.mag2();
     fTransportEndKineticEnergy  = 0.5 * restMass * veloc_sq /
                        ( 1 - veloc_sq / c_squared );   // Lorentz correction
#endif
     fTransportEndKineticEnergy = track.GetKineticEnergy();

     // fTransportEndPolarization= aFieldTrack.GetSpin(); // Not yet possible

     fParticleIsLooping = fFieldPropagator->IsParticleLooping();
     endpointDistance= (fTransportEndPosition-startPosition).mag();
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  //  then a boundary will also limit the step -> we must flag this.
  if (currentMinimumStep == 0.0 ) {
     if( currentSafety == 0.0 ){
	fGeometryLimitedStep= true;
     }
  }

  // Update the safety starting from the end-point, if it will become 
  //  negative at the end-point.
  // 
  if( currentSafety < endpointDistance ) {
      G4double endSafety;
      endSafety = fLinearNavigator->ComputeSafety( fTransportEndPosition);
      currentSafety= endSafety;
      fPreviousSftOrigin = fTransportEndPosition;
      fPreviousSafety= currentSafety; 
      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength
      currentSafety += endpointDistance;

#ifdef G4DEBUG_TRANSPORT      
      cout.precision(5);
      cout << "***Transportation::AlongStepGPIL ** " << G4endl ;
      cout << "  Called Navigator->ComputeSafety " << G4endl
	   << "    with position = " << fTransportEndPosition << G4endl
	   << "    and it returned safety= " << endSafety << G4endl; 
      cout << "  I add the endpoint distance " << endpointDistance 
	   << "   to it " 
	   << "   to obtain a pseudo-safety= " << currentSafety 
	   << "   which I return."  << G4endl; 
#endif
  }				    

  fParticleChange.SetTrueStepLength(geometryStepLength) ;

  return geometryStepLength;
}


G4VParticleChange* G4Transportation::AlongStepDoIt(
			     const G4Track& track,
			     const G4Step&  stepData
			    )
{
  //   Initialize ParticleChange  (by setting all its members equal
  //                               to corresponding members in G4Track)
  fParticleChange.Initialize(track);

  //
  //  Code for specific process 
  
  fParticleChange.SetPositionChange(fTransportEndPosition);
  fParticleChange.SetMomentumChange(fTransportEndMomentumDir);
  fParticleChange.SetEnergyChange(fTransportEndKineticEnergy);
  fParticleChange.SetMomentumChanged(fMomentumChanged);

  G4double deltaTime=0.0;
#if HARMONIC_MEAN_VELOCITY
  G4double meanInverseVelocity;
  meanInverseVelocity= 0.5/stepData.GetPreStepPoint()->GetVelocity()+
                       0.5/stepData.GetPostStepPoint()->GetVelocity();  
  if ( meanInverseVelocity < kInfinity ) {
     deltaTime= track.GetStepLength() * meanInverseVelocity; 
  }
#endif
  G4double finalVelocity= track.GetVelocity();
  if ( finalVelocity > 0.0 ) {
     deltaTime= track.GetStepLength() / finalVelocity; 
  }  

  fParticleChange. SetTimeChange( track.GetGlobalTime() + deltaTime );

  // Now Correct by Lorentz factor to get "proper" deltaTime
  //
  G4double  restMass = track.GetDynamicParticle()->GetMass();
  G4double deltaProperTime= deltaTime * (restMass / track.GetTotalEnergy());

  fParticleChange. SetProperTimeChange(track.GetProperTime() 
				                        + deltaProperTime );
  // fParticleChange.SetEnergyChange( Energy );
  //fParticleChange. SetTrueStepLength( track.GetStepLength() );

#ifdef DETECT_LOOPER
  // If the particle is caught looping in a magnetic field (doing many steps)
  //    this kills it ...
  // But currently a user-limit maximum Step size alleviates this problem,
  //    so this code is no longer used.
  if ( fParticleIsLooping ){
      // Kill the looping particle  
      fParticleChange.SetStatusChange( fStopAndKill ) ;
      // ClearNumberOfInteractionLengthLeft();
  }
#endif

  return &fParticleChange;

}

// This ensures that the PostStep action is always called,
//   so that it can do the relocation if it is needed.
// 
G4double 
G4Transportation::PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4double   previousStepSize,
			     G4ForceCondition* pForceCond
			    )
{ 
  *pForceCond= Forced; 

  return DBL_MAX;  // was kInfinity; but convention now is DBL_MAX
}

G4VParticleChange* G4Transportation::PostStepDoIt(
			     const G4Track& track,
			     const G4Step&  stepData
			    )
{
  const G4VTouchable* retCurrentTouchable;   // The one to return

  //   Initialize ParticleChange  (by setting all its members equal
  //                               to corresponding members in G4Track)
  // 
  // fParticleChange.Initialize(track);  // To initialise TouchableChange
  fParticleChange.SetStatusChange(track.GetTrackStatus());

  // If the Step was determined by the volume boundary,
  // logically relocate the particle
  //
  if( fGeometryLimitedStep ){  
    // fCurrentTouchable will now become the previous touchable, 
    //  and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)

    SetTheOtherTouchableFree(fCurrentTouchable);
    fCurrentTouchable= GetFreeTouchable();

    fLinearNavigator->SetGeometricallyLimitedStep();
    fLinearNavigator-> LocateGlobalPointAndUpdateTouchable( 
                                                track.GetPosition(),
						track.GetMomentumDirection(),
						fCurrentTouchable,
			                        true);

    // Check whether the particle is out of the world volume 
    //   If so it has exited and must be killed.
    if( fCurrentTouchable->GetVolume() == 0 ){
       fParticleChange.SetStatusChange( fStopAndKill ) ;
    }
    retCurrentTouchable= fCurrentTouchable;
    fParticleChange.SetTouchableChange( fCurrentTouchable );
  }
  else{                    // fGeometryLimitedStep  is false
#ifdef G4DEBUG
    // Although the location is changed, we know that the physical 
    //   volume remains constant. 
    // In order to help in checking the user geometry
    //    we perform a full-relocation and check its result 
    //     *except* if we have made a very small step from a boundary
    //      (ie remaining inside the tolerance

    G4bool  startAtSurface_And_MoveEpsilon;
    startAtSurface_And_MoveEpsilon=
             (stepData.GetPreStepPoint()->GetSafety() == 0.0)
          && (stepData.GetStepLength() < kCarTolerance );
    if( startAtSurface_And_MoveEpsilon) {

       // fCurrentTouchable will now become the previous touchable, 
       SetTheOtherTouchableFree(fCurrentTouchable);
       fCurrentTouchable= GetFreeTouchable();

       fLinearNavigator-> LocateGlobalPointAndUpdateTouchable( 
                                                track.GetPosition(),
						track.GetMomentumDirection(),
						fCurrentTouchable,
			                        true);
       if( fCurrentTouchable->GetVolume() != track.GetVolume() ){
          // 
          G4cerr << " ERROR: A relocation within safety has caused a volume change! " << G4endl ; 
          G4cerr << "   The old volume is called " 
	         << track.GetVolume()->GetName() << G4endl; 
          G4cerr << "   The new volume is called ";
          if ( fCurrentTouchable->GetVolume() != 0 )
	     G4cerr << fCurrentTouchable->GetVolume()->GetName() << G4endl; 
          else
	     G4cerr << "Out of World" << G4endl; 

          G4cerr.precision(7);
          G4cerr << "   The position is " << track.GetPosition() <<  G4endl;

          // Let us relocate again, for debuging
          fLinearNavigator-> LocateGlobalPointAndUpdateTouchable( 
                                                track.GetPosition(),
						track.GetMomentumDirection(),
						fCurrentTouchable,
			                        true);
          G4cerr << "   The newer volume is called " ;
          if ( fCurrentTouchable->GetVolume() != 0 )
	     G4cerr << fCurrentTouchable->GetVolume()->GetName() << G4endl; 
          else
	     G4cerr << "Out of World" << G4endl; 
       }

       assert( fCurrentTouchable->GetVolume()->GetName() == 
               track.GetVolume()->GetName() );
       retCurrentTouchable = fCurrentTouchable; 
       fParticleChange.SetTouchableChange( fCurrentTouchable );
       
    }else{
       retCurrentTouchable = track.GetTouchable();
       fParticleChange.SetTouchableChange( track.GetTouchable() );
    }
    //  This must be done in the above if ( AtSur ) fails
    //  We also do it for if (true) in order to get debug/opt to  
    //  behave as exactly the same way as possible.
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition());
#else
    // ie #ifndef G4DEBUG does a quick relocation
    
    // The serves only to move the Navigator's location
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition());

    // The value of the track's current Touchable is retained. 
    //    (and it must be correct because we must use it below to
    //      overwrite the (unset) one in particle change)
    //  Although in general this is fCurrentTouchable, at the start of
    //   a step it could be different ... ??
    fParticleChange.SetTouchableChange( track.GetTouchable() );
    retCurrentTouchable = track.GetTouchable();
#endif

  }                   // endif ( fGeometryLimitedStep ) 

  const G4VPhysicalVolume *pNewVol = retCurrentTouchable->GetVolume();
  const G4Material *pNewMaterial=0; 
  if( pNewVol != 0 ) pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial(); 

  // ( <const_cast> pNewMaterial );
  fParticleChange.SetMaterialChange( (G4Material *) pNewMaterial );
  //    temporarily until Get/Set Material of ParticleChange, 
  //    and StepPoint can be made const. 


  // Set the touchable in ParticleChange
  //   this must always be done because the particle change always
  //   uses this value to overwrite the current touchable pointer.
  //
  fParticleChange.SetTouchableChange(retCurrentTouchable);

  return &fParticleChange;
}
