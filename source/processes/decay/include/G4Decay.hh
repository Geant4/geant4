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
// $Id: G4Decay.hh 105727 2017-08-16 12:47:05Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
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
//   PreAssignedDecayTime         18 Jan. 2001 H.Kurashige
//   Add External Decayer         23 Feb. 2001  H.Kurashige
//   Remove PhysicsTable          12 Feb. 2002 H.Kurashige
//   Fixed bug in PostStepGPIL 
//    in case of stopping during AlongStepDoIt 12 Mar. 2004 H.Kurashige
//   Add GetRemainderLifeTime  10 Aug/2004 H.Kurashige
//   Add DaughterPolarization     23 July 2008 H.Kurashige


#ifndef G4Decay_h
#define G4Decay_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4DecayProcessType.hh"

class G4VExtDecayer;

class G4Decay : public G4VRestDiscreteProcess 
{
 // Class Description
  //  This class is a decay process

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

  public: //With Description
     // G4Decay Process has both 
     // PostStepDoIt (for decay in flight) 
     //   and 
     // AtRestDoIt (for decay at rest)
  
     virtual G4VParticleChange *PostStepDoIt(
			     const G4Track& aTrack,
                             const G4Step& aStep
                            ) override;

     virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    ) override;

     virtual void BuildPhysicsTable(const G4ParticleDefinition&) override; 
     // In G4Decay, thePhysicsTable stores values of
    //    beta * std::sqrt( 1 - beta*beta) 
    //  as a function of normalized kinetic enregy (=Ekin/mass),
    //  becasuse this table is universal for all particle types,


    virtual G4bool IsApplicable(const G4ParticleDefinition&) override;
    // returns "true" if the decay process can be applied to
    // the particle type. 
 
  protected: // With Description
    virtual G4VParticleChange* DecayIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    );
    // The DecayIt() method returns by pointer a particle-change object,
    // which has information of daughter particles.

    // Set daughter polarization
    //  NO OPERATION in the base class of G4Decay 
    virtual void DaughterPolarization(const G4Track& aTrack,
			      G4DecayProducts* products);

 public:
    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            ) override;

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            ) override;

  protected: // With Description
    // GetMeanFreePath returns ctau*beta*gamma for decay in flight 
    // GetMeanLifeTime returns ctau for decay at rest
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double   previousStepSize,
                              G4ForceCondition* condition
                             ) override;

    virtual G4double GetMeanLifeTime(const G4Track& aTrack,
                              G4ForceCondition* condition
                            ) override;

   public: //With Description
     virtual void StartTracking(G4Track*) override;
     virtual void EndTracking() override;
      // inform Start/End of tracking for each track to the physics process 

   public: //With Description
     void SetExtDecayer(G4VExtDecayer*);
     const G4VExtDecayer* GetExtDecayer() const;
     // Set/Get External Decayer
   
    G4double GetRemainderLifeTime() const;  
    //Get Remainder of life time at rest decay 

    virtual void ProcessDescription(std::ostream& outFile) const override;
    //

  protected:
     G4int verboseLevel;
     // controle flag for output message
     //  0: Silent
     //  1: Warning message
     //  2: More

  protected:
    // HighestValue.
    const G4double HighestValue;
 
    // Remainder of life time at rest
    G4double                 fRemainderLifeTime;
  
    // ParticleChange for decay process
    G4ParticleChangeForDecay fParticleChangeForDecay;
    
    // External Decayer
    G4VExtDecayer*    pExtDecayer;
};

inline  
  G4VParticleChange* G4Decay::AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}


inline
 const G4VExtDecayer* G4Decay::GetExtDecayer() const
{
  return pExtDecayer;
}

inline
 G4double G4Decay::GetRemainderLifeTime() const 
{
  return fRemainderLifeTime;
}

#endif










