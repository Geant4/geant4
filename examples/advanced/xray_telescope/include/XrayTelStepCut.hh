// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelStepCut.hh                               *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
//
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelStepCut_h
#define XrayTelStepCut_h 1

#include "G4ios.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Step.hh"
#include "globals.hh"

class XrayTelStepCut : public G4VDiscreteProcess
{
public:     

  XrayTelStepCut(const G4String& processName ="UserStepCut" );
  XrayTelStepCut(XrayTelStepCut &);

  ~XrayTelStepCut();

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
  XrayTelStepCut & operator=(const XrayTelStepCut &right);

private:

  G4double MaxChargedStep ;
};

// inlined function members implementation

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

inline G4double XrayTelStepCut::PostStepGetPhysicalInteractionLength(
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

inline G4VParticleChange* XrayTelStepCut::PostStepDoIt(
						       const G4Track& aTrack,
						       const G4Step&
						       )
{
  // do nothing
  aParticleChange.Initialize(aTrack);
  return &aParticleChange;
}

inline G4double XrayTelStepCut::GetMeanFreePath(const G4Track& aTrack,
						G4double   previousStepSize,
						G4ForceCondition* condition
						)
{
  return 0.;
}

#endif

