// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Transportation.cc,v 1.2 1999-01-08 11:23:58 gunter Exp $
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
// Created:  19 March 1997, J. Apostolakis
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
  G4double geometryStepLength, endpointSafety=0;

  fParticleIsLooping = false;

  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value
  *selection = CandidateForSelection;

  // Get initial Energy/Momentum of the track
  //
  const G4DynamicParticle*  pParticle  = track.GetDynamicParticle();
  // const G4double      startEnergy      = pParticle->GetKineticEnergy();
  const G4ThreeVector startMomentumDir = pParticle->GetMomentumDirection();
  const G4ThreeVector startPosition  = track.GetPosition();
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
        G4double newSafety;

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

     endpointSafety= currentSafety - endpointDistance;
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

     // Remember last safety origin & value.
     fPreviousSftOrigin = startPosition;
     fPreviousSafety= currentSafety; 		    
        
     // Get the End-Position and End-Momentum (Dir-ection)
     fTransportEndPosition= aFieldTrack.GetPosition();
     // Momentum:  Magnitude and direction can be changed too now ...
     fMomentumChanged = true; 
     fTransportEndMomentumDir= aFieldTrack.GetMomentumDir();

     // fTransportEndKineticEnergy= aFieldTrack.GetEnergy(); // Energy is wrong
#if EBUG_ENERGY_IN_FIELD
     G4ThreeVector endVelocity = aFieldTrack.GetVelocity();
     G4double  veloc_sq = endVelocity.mag2();
     fTransportEndKineticEnergy  = 0.5 * restMass * veloc_sq /
                       ( 1. - (veloc_sq / c_squared) );  // Lorentz correction
#endif
     fTransportEndKineticEnergy = track.GetKineticEnergy();

     // fTransportEndPolarization= aFieldTrack.GetSpin(); // Not yet possible

     fParticleIsLooping = fFieldPropagator->IsParticleLooping();
     endpointDistance= (fTransportEndPosition-startPosition).mag();

     // Recover an endpoint safety
     fPreviousSftOrigin= fFieldPropagator->GetLastSafetyOrigin();
     fPreviousSafety=    fFieldPropagator->GetLastSafetyValue();

     // Calculate the endpoint safety
     G4double distEnd2LastSafOrigin;
     distEnd2LastSafOrigin= (fTransportEndPosition-fPreviousSftOrigin).mag();
     endpointSafety = fPreviousSafety - distEnd2LastSafOrigin;
  }

  // If we are asked to go a step length of 0, and we are on a boundary
  //  then a boundary will also limit the step -> we must flag this.
  if (currentMinimumStep == 0.0 ) {
     if( currentSafety == 0.0 ){
	fGeometryLimitedStep= true;
     }
  }
  // Check the endpointSafety
  G4double distEnd2Start;
  distEnd2Start= (fTransportEndPosition-startPosition).mag();
  endpointSafety = fPreviousSafety - distEnd2Start;

  // Update the safety starting from the end-point, if it will become 
  //  negative at the end-point.
  // 
  if( endpointSafety < 0 ){
      endpointSafety = fLinearNavigator->ComputeSafety( fTransportEndPosition);
      fPreviousSftOrigin = fTransportEndPosition;
      fPreviousSafety= endpointSafety; 

#ifdef G4DEBUG_TRANSPORT      
      cout.precision(5);
      cout << "***Transportation::AlongStepGPIL ** " << endl ;
      cout << "  Called Navigator->ComputeSafety " << endl
	   << "    with position = " << fTransportEndPosition << endl
	   << "    and it returned safety= " << endpointSafety << endl; 
#endif
  }				    

  // Because the Stepping Manager assumes it is from the start point, 
  //  we must add the StepLength to create a "pseudo" safety at start point.
  currentSafety = endpointDistance + endpointSafety;

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
  G4VTouchable* retCurrentTouchable;   // The one to return

  //   Initialize ParticleChange  (by setting all its members equal
  //                               to corresponding members in G4Track)
  // 
  // fParticleChange.Initialize(track);  // To initialise TouchableChange
  fParticleChange.SetStatusChange(track.GetTrackStatus());

  // fCurrentTouchable will now become the previous touchable, 
  //  and what was the previous will be freed.
  // (We need this because the preStepPoint can point to the previous 
  //  touchable)
  //
  // SetTheOtherTouchableFree(fCurrentTouchable);  // Do it only if needed.
  // fCurrentTouchable= GetFreeTouchable();        // Do it only if needed.

  // If the Step was determined by the volume boundary,
  // logically relocate the particle
  // if( stepData.GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
  // If the use of fGeomBoundary is suppressed, we can probably change this to:
  //  if( stepData.GetPostStepPoint()->GetProcessDefinedStep() == this ){
  //
  if( fGeometryLimitedStep ){  
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
  }
  else{
#ifdef G4VERBOSE
    // fCurrentTouchable will now become the previous touchable, 
    SetTheOtherTouchableFree(fCurrentTouchable);
    fCurrentTouchable= GetFreeTouchable();

    // Although the location is changed, we know that the physical 
    //   volume remains constant. 
    // Currently a pseudo-relocation is/was required here:
    fLinearNavigator-> LocateGlobalPointAndUpdateTouchable( 
                                                track.GetPosition(),
						track.GetMomentumDirection(),
						fCurrentTouchable,
			                        true);
    if( fCurrentTouchable->GetVolume() != track.GetVolume() ){
       // 
       G4cerr << " ERROR: A relocation within safety has caused a volume change! " << endl ; 
       G4cerr << "   The old volume is called " 
	      << track.GetVolume()->GetName() << endl; 
       G4cerr << "   The new volume is called ";
       if ( fCurrentTouchable->GetVolume() != 0 )
	  G4cerr << fCurrentTouchable->GetVolume()->GetName() << endl; 
       else
	  G4cerr << "Out of World" << endl; 

       G4cerr.precision(7);
       G4cerr << "   The position is " << track.GetPosition() <<  endl;

       // Let us relocate again, for debuging
       fLinearNavigator-> LocateGlobalPointAndUpdateTouchable( 
                                                track.GetPosition(),
						track.GetMomentumDirection(),
						fCurrentTouchable,
			                        true);
       G4cerr << "   The newer volume is called " ;
       if ( fCurrentTouchable->GetVolume() != 0 )
	  G4cerr << fCurrentTouchable->GetVolume()->GetName() << endl; 
       else
	  G4cerr << "Out of World" << endl; 
    }

    assert( fCurrentTouchable->GetVolume()->GetName() == 
            track.GetVolume()->GetName() );
    retCurrentTouchable = fCurrentTouchable; 
#else
    // ie #ifndef G4VERBOSE does a quick relocation
    
    // The serves only to move the Navigator's location
    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition());
    // The value of the track's current Touchable is retained. 
    //    (and it must be correct because we must use it below to
    //      overwrite the (unset) one in particle change)
    //  Although in general this is fCurrentTouchable, at the start of
    //   a step it could be different ... ??
    retCurrentTouchable = track.GetTouchable();
#endif

  }

  // Set the touchable in ParticleChange
  //   this must always be done because the particle change always
  //   uses this value to overwrite the current touchable pointer.
  //
  fParticleChange.SetTouchableChange(retCurrentTouchable);

  return &fParticleChange;
}
