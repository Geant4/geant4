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
/// \file field/field03/include/F03StepCut.hh
/// \brief Definition of the F03StepCut class
//
// $Id$
// 

#ifndef F03StepCut_h
#define F03StepCut_h 1

#include "G4VDiscreteProcess.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "globals.hh"

class F03StepCut : public G4VDiscreteProcess
{
  public:     

     F03StepCut(const G4String& processName ="UserStepCut");
     F03StepCut(F03StepCut &);

     ~F03StepCut();

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

     virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );

     void SetMaxStep(G4double);

  protected:

     // it is not needed here !
     virtual G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

                            
  private:
  
     // hide assignment operator as private 
     F03StepCut & operator=(const F03StepCut &right);

  private:

     G4double fMaxChargedStep ;
};

// inlined function members implementation

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

inline G4double F03StepCut::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
                             G4double,
                             G4ForceCondition* condition
                            )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double proposedStep = DBL_MAX;

   if((fMaxChargedStep > 0.) &&
      (aTrack.GetVolume() != 0) &&
      (aTrack.GetVolume()->GetName() == "Absorber") &&
      (aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge() != 0.))
        proposedStep = fMaxChargedStep ;

   return proposedStep;
}

inline G4VParticleChange* F03StepCut::PostStepDoIt(
                             const G4Track& aTrack,
                             const G4Step&
                            )
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

inline G4double F03StepCut::GetMeanFreePath(const G4Track&,
                             G4double,
                             G4ForceCondition*)
{
  return 0.;
}

#endif
