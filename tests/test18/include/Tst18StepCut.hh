// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18StepCut.hh,v 1.4 2000-06-14 12:33:13 flei Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// ------------------------------------------------------------
//                  05 March 1999  L.Urban
// ------------------------------------------------------------

#ifndef Tst18StepCut_h
#define Tst18StepCut_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Step.hh"

class Tst18StepCut : public G4VDiscreteProcess
{
  public:     

     Tst18StepCut(const G4String& processName ="UserStepCut" );
     Tst18StepCut(Tst18StepCut &);

     ~Tst18StepCut();

     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

    void SetMaxStep(G4double);

  protected:

     // it is not needed here !
     G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

			    
  private:
  
  // hide assignment operator as private 
      Tst18StepCut & operator=(const Tst18StepCut &right);

  private:

     G4double MaxChargedStep ;
};

// inlined function members implementation

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

inline G4double Tst18StepCut::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double ProposedStep = DBL_MAX;

   if((MaxChargedStep > 0.) &&
      (aTrack.GetVolume() != NULL) &&
      (aTrack.GetVolume()->GetName() == "Absorber") &&
      (aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge() != 0.))
        ProposedStep = MaxChargedStep ;

   return ProposedStep;
}

inline G4VParticleChange* Tst18StepCut::PostStepDoIt(
                             const G4Track& aTrack,
                             const G4Step&
                            )
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

inline G4double Tst18StepCut::GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
  return 0.;
}

#endif

