// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransportationIntegralTOF.hh,v 1.2 1999-12-15 14:53:50 gunter Exp $
// Started from  G4Transportation.hh,v 2.4 1998/11/11 20:03:58 japost 
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
#ifndef G4TransportationIntegralTOF_hh
#define G4TransportationIntegralTOF_hh 1

#include "G4VProcess.hh"

class G4Navigator;
class G4PropagatorInField;
class G4Track;
class G4Step;
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleChangeForTransport.hh"

class G4TransportationIntegralTOF : public G4VProcess 
{
  // Concrete class that does the geometrical transport 

  public:
     G4TransportationIntegralTOF();
     ~G4TransportationIntegralTOF(); 

     //  G4double          GetContinuousStepLimit  (
     G4double      AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
				   G4double  previousStepSize,
			           G4double  currentMinimumStep, 
				   G4double& currentSafety,
				   G4GPILSelection* selection
			    );

     G4VParticleChange* AlongStepDoIt(
			     const G4Track& track,
			     const G4Step& stepData
			    );

     // This only does the relocation
     // 
     G4VParticleChange* PostStepDoIt(
			     const G4Track& track,
			     const G4Step&  stepData
			    );

     // This forces the PostStepDoIt action to be called, 
     //   but does not limit the step.
     //
     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4double   previousStepSize,
			     G4ForceCondition* pForceCond
			    );

     //  Access/set the assistant class that Propagate in a Field
     G4PropagatorInField* GetPropagatorInField();
     void SetPropagatorInField( G4PropagatorInField* pFieldPropagator);

     //  no operation in  AtRestDoIt
     G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  AtRestDoIt
     G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return NULL;};
  protected:
     //  Checks whether a field exists for the "global" field manager.
     G4bool               DoesGlobalFieldExist();

  private:
     // Needed to add the Relocation in the PostStepDoIt
     //                                         Copied from SteppingManager
     G4VTouchable* GetFreeTouchable();
        // Get Touchable which is free, i.e. not assigined to Track/StepPoint
        // If no free touchable is availabe, the NULL will be returned
        // Once you get a Touchable, it will be set to NotFree. 
     void SetTheOtherTouchableFree(G4VTouchable* pTouch);
        // Set the partner of the given Touchable to be free. For example,
        // the argument has fTouchable1, then fTouchable2 will be set to 

  private:
     // The Propagators used to transport the particle
     G4Navigator*         fLinearNavigator;
     G4PropagatorInField* fFieldPropagator;

     // Field Manager for the whole Detector
     // G4FieldManager*      fGlobalFieldMgr;     // Used MagneticField CC

     // The particle's state after this Step, Store for DoIt
     G4ThreeVector        fTransportEndPosition;
     G4ThreeVector        fTransportEndMomentumDir;
     G4double             fTransportEndKineticEnergy;
     G4ThreeVector        fTransportEndSpin;
     G4bool               fMomentumChanged;
     G4bool               fEnergyChanged;

     G4bool               fParticleIsLooping;

     G4VTouchable*        fCurrentTouchable;
     
     // Whether a magnetic field exists ...
     // G4bool         fFieldExists;
     //   The above data member is problematic: it is useful only if
     // it is initialised.  However the transportation process(es) are not
     // capable of doing this initialisation itself (themselves) and there
     // seems no alternative agent capable of doing it right now.
     // Eg, at construction time it is likely that the field manager has
     // not yet been informed about the detector's field 
     // I cannot foresee how the transportation can be informed later. JA 
     //   The current answer is to ignore this data member and use 
     // the member function DoesGlobalFieldExist() in its place ...
     //    John Apostolakis, July 7, 1997

     // Needed for the relocation - Copied from SteppingManager
     G4VTouchable* fTouchable1;
     G4VTouchable* fTouchable2;
     G4bool fIsTouchable1Free;
     G4bool fIsTouchable2Free;

     // Flag to determine whether a boundary was reached.
     G4bool fGeometryLimitedStep;

     // Remember last safety origin & value.
     G4ThreeVector  fPreviousSftOrigin;
     G4double       fPreviousSafety; 

     // New ParticleChange
     G4ParticleChangeForTransport fParticleChange;

     G4double endpointDistance;
};


#include "G4TransportationIntegralTOF.icc"

#endif  

// End of G4TransportationIntegralTOF.hh
