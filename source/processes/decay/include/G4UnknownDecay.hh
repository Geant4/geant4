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
// $Id: G4UnknownDecay.hh 105727 2017-08-16 12:47:05Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
#ifndef G4UnknownDecay_h
#define G4UnknownDecay_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"

class G4UnknownDecay : public G4VDiscreteProcess 
{
 // Class Description
  //  This class is a decay process for "unknown" particle.

  public:
    //  Constructors 
    G4UnknownDecay(const G4String& processName ="UnknownDecay");

    //  Destructor
    virtual ~G4UnknownDecay();

  private:
    //  copy constructor
      G4UnknownDecay(const G4UnknownDecay &right);

    //  Assignment Operation (generated)
      G4UnknownDecay & operator=(const G4UnknownDecay &right);

  public: //With Description
  
     virtual G4VParticleChange *PostStepDoIt(
			     const G4Track& aTrack,
                             const G4Step& aStep
                            ) override;

     virtual void BuildPhysicsTable(const G4ParticleDefinition&) override; 
     // In G4UnknownDecay, thePhysicsTable stores values of
    //    beta * std::sqrt( 1 - beta*beta) 
    //  as a function of normalized kinetic enregy (=Ekin/mass),
    //  becasuse this table is universal for all particle types,


    virtual G4bool IsApplicable(const G4ParticleDefinition&) override;
    // returns "true" if the decay process can be applied to
    // the particle type. 
 
    virtual void ProcessDescription(std::ostream& outFile) const override;
    //

  protected: // With Description
    virtual G4VParticleChange* DecayIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    ) ;
    // The DecayIt() method returns by pointer a particle-change object,
    // which has information of daughter particles.

  public:

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            ) override;


  protected: // With Description
    // GetMeanFreePath returns ctau*beta*gamma for decay in flight 
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double   previousStepSize,
                              G4ForceCondition* condition
                             ) override;

  private:
     G4int verboseLevel;
     // controle flag for output message
     //  0: Silent
     //  1: Warning message
     //  2: More

  private:
    // HighestValue.
    const G4double HighestValue;
 
    // ParticleChange for decay process
    G4ParticleChangeForDecay fParticleChangeForDecay;
    
};

inline 
 G4double G4UnknownDecay::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   /*previousStepSize*/,
                             G4ForceCondition* condition
                            )
{
  // pre-assigned UnknownDecay time
  G4double pTime = track.GetDynamicParticle()->GetPreAssignedDecayProperTime();

  if (pTime < 0.) pTime = DBL_MIN;

  // condition is set to "Not Forced"
  *condition = NotForced;
  
  // reminder proper time
  G4double remainder = pTime - track.GetProperTime();
  if (remainder <= 0.0) remainder = DBL_MIN;
  
  // use pre-assigned Decay time to determine PIL
  //return GetMeanFreePath(track, previousStepSize, condition);
  return remainder*CLHEP::c_light;

}

inline  
  G4VParticleChange* G4UnknownDecay::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}

#endif










