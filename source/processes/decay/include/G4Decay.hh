// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Decay.hh,v 1.3 1999-12-15 14:51:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      7 July 1996 H.Kurashige
// ------------------------------------------------------------
//  New Physics scheme           18 Jan. 1997  H.Kurahige
// ------------------------------------------------------------
//   modified                     4  Feb. 1997  H.Kurahige
//   modified                     8  Sep. 1997  H.Kurahige
//   remove BuildPhysicsTable()   27 Nov. 1997   H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige
//   added aPhysicsTable          2  Aug. 1998 H.Kurashige

#ifndef G4Decay_h
#define G4Decay_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"

class G4Decay : public G4VRestDiscreteProcess 
{
  public:
    //  Constructors 
    G4Decay(const G4String& processName ="Decay");

    //  Destructor
    virtual ~G4Decay();

  private:
    //  copy constructor
      G4Decay(const G4Decay &right);

    //  Assignment Operation (generated)
      G4Decay & operator=(const G4Decay &right);

  public:
     // G4Decay Process has both 
     // PostStepDoIt (for decay in flight) 
     //   and 
     // AtRestDoIt (for decay at rest)
  
     virtual G4VParticleChange *PostStepDoIt(
			     const G4Track& aTrack,
                             const G4Step& aStep
                            );

     virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    );

     virtual void BuildPhysicsTable(const G4ParticleDefinition&); 
     // In G4Decay, thePhysicsTable stores values of
    //    beta * sqrt( 1 - beta*beta) 
    //  as a function of normalized kinetic enregy (=Ekin/mass),
    //  becasuse this table is universal for all particle types,


    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    // returns "true" if the decay process can be applied to
    // the particle type. 
 
  protected:
    virtual G4VParticleChange* DecayIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    );
    // The DecayIt() method returns by pointer a particle-change object,
    // which has information of daughter particles.

  public:
    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            );

  protected:
    // GetMeanFreePath returns ctau*beta*gamma for decay in flight 
    // GetMeanLifeTime returns ctau for decay at rest
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double   previousStepSize,
                              G4ForceCondition* condition
                             );

    virtual G4double GetMeanLifeTime(const G4Track& aTrack,
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
    // In G4Decay, thePhysicsTable stores values of
    //    beta * sqrt( 1 - beta*beta) 
    //  as a function of normalized kinetic enregy (=Ekin/mass)
    // TotBin = number of bins in thePhysicsTable
    // The PhysicsTable is created for the range of  
    // from LowestBinValue to HighestBinValue.
    const G4double LowestBinValue;
    const G4double HighestBinValue;
    const G4int    TotBin;
    
    G4PhysicsTable* aPhysicsTable;
 
    // Remainder of life time at rest
    G4double                 fRemainderLifeTime;
 
    // ParticleChange for decay process
    G4ParticleChangeForDecay fParticleChangeForDecay;
};

inline
  G4double G4Decay::AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            )
{
  fRemainderLifeTime = 
    G4VRestDiscreteProcess::AtRestGetPhysicalInteractionLength(
                             track, condition );
  return fRemainderLifeTime;
}

inline
 void  G4Decay::SetVerboseLevel(G4int value){ verboseLevel = value; }

inline
 G4int G4Decay::GetVerboseLevel() const { return verboseLevel; }

inline  
  G4VParticleChange* G4Decay::AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}

inline  
  G4VParticleChange* G4Decay::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}


#endif










