// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Started from G4Transportation.cc,v 2.18 1998/12/14 18:27:32 japost 
//
// $ Id:  $
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
//   This class' object is a process responsible for the transportation of 
// a particle, using an integral approach for estimating the Time of Flight
// for a charged particle.  Once it work, the will be incorporated again
// into G4Transportation, and this will be suppressed (deleted.)
//
// =======================================================================
// Created:  17 December 1998, J. Apostolakis
// =======================================================================

#include "G4Transportation.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4EnergyLossTables.hh"

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
      cout << "***Transportation::AlongStepGPIL ** " << endl ;
      cout << "  Called Navigator->ComputeSafety " << endl
	   << "    with position = " << fTransportEndPosition << endl
	   << "    and it returned safety= " << endSafety << endl; 
      cout << "  I add the endpoint distance " << endpointDistance 
	   << "   to it " 
	   << "   to obtain a pseudo-safety= " << currentSafety 
	   << "   which I return."  << endl; 
#endif
  }				    

  return geometryStepLength;
}

// ----------------------------------------------------------------------------
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

  // fParticleChange.SetEnergyChange( track.GetKineticEnergy() );
  //fParticleChange. SetTrueStepLength( track.GetStepLength() );

  G4double TruePathLength= track.GetStepLength();
  G4double deltaTime=0.0;

  G4double finalVelocity= track.GetVelocity();
  if ( finalVelocity > 0.0 ) {
     deltaTime= TruePathLength / finalVelocity; 
  }  

#if HARMONIC_MEAN_VELOCITY
  G4double meanInverseVelocity;
  meanInverseVelocity= 0.5/stepData.GetPreStepPoint()->GetVelocity()+
                       0.5/stepData.GetPostStepPoint()->GetVelocity();  
  if ( meanInverseVelocity < kInfinity ) {
     deltaTime= track.GetStepLength() * meanInverseVelocity; 
  }
#endif

  fParticleChange. SetTimeChange( track.GetGlobalTime() + deltaTime );
  if( 1 ) {  
  // if( verbose > 2 ) {
     G4double LabTime_mean=0.0;
     G4double velocity_start, velocity_end;  // Start/End of Step
     
     velocity_start= stepData.GetPreStepPoint() ->GetVelocity();
     velocity_end  = stepData.GetPostStepPoint()->GetVelocity(); 

     // Harmonic mean of velocities:
     //   1/harm_mean = 0.5 * ( 1/v1 + 1/v2 )
     //
     G4double meanInverseVelocity;
     if( velocity_end != 0 ) {
        meanInverseVelocity= 0.5 * (velocity_start + velocity_end) /
                                   (velocity_start * velocity_end);
        LabTime_mean= TruePathLength * meanInverseVelocity; 
     }else{
        //  The particle is comming to rest.
        //  So this estimate is completely inaccurate 
        //                                      - but how can we do better?
        G4double meanVelocity= 0.5*velocity_start;   
        LabTime_mean= TruePathLength / meanVelocity; 
     }

     if (LabTime_mean!=0.0) {
        cout << setw(12) << (deltaTime/LabTime_mean - 1.0);
     }
     cout << endl;
  }     
  // Now Correct by Lorentz factor to get "proper" deltaTime
  //
  G4double  restMass = track.GetDynamicParticle()->GetMass();
  G4double deltaProperTime= deltaTime * (restMass / track.GetTotalEnergy());

  fParticleChange. SetProperTimeChange(track.GetProperTime() 
		 		                          + deltaProperTime );

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

// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------

G4VParticleChange* G4Transportation::PostStepDoIt(
			     const G4Track& track,
			     const G4Step&  stepData
			    )
{
  G4VTouchable* retCurrentTouchable;   // The one to return

  //   Initialize ParticleChange  (by setting all its members equal
  //                               to corresponding members in G4Track)
  fParticleChange.Initialize(track); //  Initialises Times
                                     //   and        TouchableChange
  // fParticleChange.SetStatusChange(track.GetTrackStatus()); // If not above Init
                                                              // MUST do this anyway

  static G4int fVerbose = 0;  // 0=nothing,  1-5=tabular,  6=very verbose
#ifndef NO_CHECK
  // check
  G4double DeltaPosition= (fTransportEndPosition - stepData.GetPreStepPoint()
			                                   ->GetPosition()).mag();
  static G4bool checkDeltaPos=0; // To switch on/off the check when debugging
  cout.precision(8);

  if ( fVerbose > 5 ) {
    cout << "G4Transportation::PostStepDoIt ***  " << endl
	 << " The track has gone a geometrical distance of " << endl
	 << "   new estimate=" <<  DeltaPosition << "  (PostStepDoIt) " << endl
	 << "   previously  =" <<  endpointDistance 
	 <<    "  was the distance of endpoint (in AlongStepGPIL) "  << endl;
    cout.precision(4);
    cout << "   difference    =" <<  DeltaPosition-endpointDistance << endl;
    if( endpointDistance != 0.0 )
      cout << "   relative diff =" <<  DeltaPosition/endpointDistance - 1.0 
	   << endl;
  }else if ( fVerbose > 0 ) {
      cout.precision(3);
      cout << setw(12) << "DeltaPosition" << " "     
	   << setw(16) << "EndpointDistance"  << " "     
	   << setw(12) << "Difference"    << " " 
	   << setw(12) << "RelativeDiff"  << " ";
      cout << endl;

      cout << setw(12) << DeltaPosition << " "     
	   << setw(15) << endpointDistance << " "     
	   << setw(12) << DeltaPosition-endpointDistance  << " " ;
      if( endpointDistance != 0.0 )
	cout << setw(12) << DeltaPosition/endpointDistance - 1.0;
      cout << endl;

  }

  if( checkDeltaPos ){
    // assert( DeltaPosition == endpointDistance );

    assert( fabs(DeltaPosition-endpointDistance) < kCarTolerance );
  }
#endif

  const G4int lwide=12;           // Width of printout column
  const G4int l2wide=2*lwide+2;

  const G4DynamicParticle*    pParticle    = track.GetDynamicParticle();
  const G4ParticleDefinition* pPartDef = pParticle->GetDefinition();
  G4double              particleCharge= pPartDef->GetPDGCharge(); 
  
  G4double LabTime_diff, ProperTime_diff;
  if( (particleCharge!=0.0) )
  {
     // Electrons & positrons have their proper time calculated from a table
     //
     G4double LabTime_start, LabTime_end; // , LabTime_diff;
     G4Material * pMaterial= track.GetMaterial();

     // Estimated time from start of step until all Energy is lost
     LabTime_start= 
         G4EnergyLossTables::GetLabTime( pPartDef, 
					 stepData.GetPreStepPoint()
						 ->GetKineticEnergy(), 
					 pMaterial);

     // Estimated time from end   of step until all Energy is lost
     LabTime_end=   
         G4EnergyLossTables::GetLabTime( pPartDef, 
					 stepData.GetPostStepPoint()
					         ->GetKineticEnergy(),  
					 pMaterial);

     LabTime_diff= LabTime_start - LabTime_end;
     // fParticleChange. SetTimeChange( track.GetGlobalTime() + LabTime_diff );


#ifdef G4VERBOSE
     if( fVerbose > 5 ){
       cout.precision(6);
       cout << "Transportation::PostStepDoIt reports: " << endl;
       cout << "  Energies " 
	    << "   start  = " << stepData.GetPreStepPoint()->GetKineticEnergy() 
	    << endl
	    << "     end  = " << stepData.GetPostStepPoint()->GetKineticEnergy()
	    << endl;
       cout << "  LabTime: " << endl
	    << "    start = " << LabTime_start << endl
	    << "    - end = " << LabTime_end   << endl
	    << "    diff  = " << LabTime_diff  << endl;
     }else if (fVerbose > 1 ){
       cout.precision(6);
       cout << "Transportation::PostStepDoIt reports: " << endl;
       cout << setw(l2wide) << "       Energies " 
	    << setw(l2wide+lwide+1) << "               LabTime " 
	    << endl ;
       cout << setw(lwide) << "  Start "  << " " << setw(lwide) << "  End "    << " "  
            << setw(lwide) << "  Start "  << " " << setw(lwide) << "  End "    << " "
            << setw(lwide) << "  Diff "  
            << endl;
       cout << setw(lwide) << stepData.GetPreStepPoint()->GetKineticEnergy() 
	    << setw(lwide) << stepData.GetPostStepPoint()->GetKineticEnergy()
	    << setw(lwide) << LabTime_start 
	    << setw(lwide) << LabTime_end   
	    << setw(lwide) << LabTime_diff  
 	    << endl;
    }
#endif

     G4double ProperTime_start, ProperTime_end; // , ProperTime_diff;

     // Estimated proper time from start of step until all Energy is lost
     ProperTime_start= 
         G4EnergyLossTables::GetProperTime( pPartDef, 
					    stepData.GetPreStepPoint()
						    ->GetKineticEnergy(), 
					    pMaterial);

     // Estimated proper time from end   of step until all Energy is lost
     ProperTime_end=   
         G4EnergyLossTables::GetProperTime( pPartDef, 
					    stepData.GetPostStepPoint()
					            ->GetKineticEnergy(),  
					    pMaterial);

     ProperTime_diff= ProperTime_start - ProperTime_end;
     // fParticleChange. SetProperTimeChange( track.GetProperTime()
     // 					   + ProperTime_diff      );

#ifdef G4VERBOSE
     if( fVerbose > 5 ){
       cout.precision(6);     cout << "  ProperTime: " << endl
	  << "    start = " << ProperTime_start << endl
	  << "    - end = " << ProperTime_end   << endl
	  << "    diff  = " << ProperTime_diff  << endl;
     }else if (fVerbose > 1 ){
       cout.precision(6);

       cout.precision(6);
       cout << setw(lwide) << "  ProperTime: ---------------- " << endl; 
       cout << setw(lwide) << "    Start  " 
	    << setw(lwide) << "     End  " 
	    << setw(lwide) << "     Diff   " << endl;
       cout << setw(lwide) << ProperTime_start << " "
	    << setw(lwide) << ProperTime_end   << " "
	    << setw(lwide) << ProperTime_diff  << " " 
	    << endl;
      }
#endif
  // FOR TEST ONLY  - start
  }
  if(0){
  // FOR TEST ONLY  - end

  }else{
     G4double Velocity= track.GetVelocity();
     G4double deltaTime=0.0;
     if ( Velocity > 0.0 ) {
        deltaTime= track.GetStepLength() / Velocity; 
     }

     // fParticleChange. SetTimeChange( track.GetGlobalTime() + deltaTime );

#ifdef G4VERBOSE
     G4double  restMass = track.GetDynamicParticle()->GetMass();
     G4double deltaProperTime= deltaTime * (restMass / track.GetTotalEnergy());
     if( fVerbose > 5 ){
       cout <<   "  DeltaLabTime = " << deltaTime << endl
	    <<   "  New method   = " << LabTime_diff << endl
	    <<   "   change      = " << deltaTime - LabTime_diff << endl;
       if (deltaTime!=0.0) 
	 cout << "  % change     = " << (LabTime_diff/deltaTime - 1.0) << endl;
       
       // Now the time corrected by Lorentz factor to get "proper" deltaTime
       //
       cout << "  DeltaProperTime" << endl;
       
       cout <<   "  Old method  = " << deltaProperTime << endl
	    <<   "  New method  = " << ProperTime_diff << endl
	    <<   "   change     = " << deltaProperTime - ProperTime_diff << endl;
       if (deltaProperTime!=0.0) 
	 cout << "  % change    = " << (ProperTime_diff/deltaProperTime - 1.0) << endl;
     }else if (fVerbose > 1 ){

       cout << "    DeltaLabTime " << endl;
       cout << setw(lwide) << " Old method " << " "
	    << setw(lwide) << " New method " << " "
	    << setw(lwide) << " Difference " << " " << endl;

       cout << setw(lwide) <<  deltaTime 
	    << setw(lwide) <<  LabTime_diff 
	    << setw(lwide) <<  deltaTime - LabTime_diff;
       if (deltaTime!=0.0) 
	 cout << "  % change     = " << (LabTime_diff/deltaTime - 1.0); 
       cout << endl;
       
       cout << "  DeltaProperTime " << endl;
       cout << setw(lwide) << " Old method " << " "
	    << setw(lwide) << " New method " << " "
	    << setw(lwide) << " Difference " << " " 
	    << setw(lwide) << " Relative%  " << " " 
	    << endl;
       cout << setw(lwide) << deltaProperTime << " "
	    << setw(lwide) << ProperTime_diff << " "
	    << setw(lwide) << deltaProperTime - ProperTime_diff << " ";
       if (deltaProperTime!=0.0) 
	 cout << "  % change    = " << (ProperTime_diff/deltaProperTime - 1.0);
       cout << endl;
     }
#endif       

  }  // end of if( charge != 0.0 )

  // ---------------------------------------------------------------------
  //  Do the RELOCATION of the particle:
  //     Inform the tracking's Navigator of the particle's new location
  // ---------------------------------------------------------------------
  //  This method has taken over the responsibility for the track's relocation
  //    now that it has been removed from G4SteppingManager.

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
