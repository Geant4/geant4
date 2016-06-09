//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4UnknownDecay.hh,v 1.2 2004/12/02 07:06:30 kurasige Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
#ifndef G4UnknownDecay_h
#define G4UnknownDecay_h 1

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
                            );

     virtual void BuildPhysicsTable(const G4ParticleDefinition&); 
     // In G4UnknownDecay, thePhysicsTable stores values of
    //    beta * std::sqrt( 1 - beta*beta) 
    //  as a function of normalized kinetic enregy (=Ekin/mass),
    //  becasuse this table is universal for all particle types,


    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    // returns "true" if the decay process can be applied to
    // the particle type. 
 
  protected: // With Description
    virtual G4VParticleChange* DecayIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    );
    // The DecayIt() method returns by pointer a particle-change object,
    // which has information of daughter particles.

  public:

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );


  protected: // With Description
    // GetMeanFreePath returns ctau*beta*gamma for decay in flight 
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double   previousStepSize,
                              G4ForceCondition* condition
                             );

  public:
     void  SetVerboseLevel(G4int value);
     G4int GetVerboseLevel() const;

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
  return remainder*c_light;

}

inline
 void  G4UnknownDecay::SetVerboseLevel(G4int value){ verboseLevel = value; }

inline
 G4int G4UnknownDecay::GetVerboseLevel() const { return verboseLevel; }

inline  
  G4VParticleChange* G4UnknownDecay::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}

#endif










