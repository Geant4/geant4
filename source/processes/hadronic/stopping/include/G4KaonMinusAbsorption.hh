// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KaonMinusAbsorption.hh,v 1.3 2000-12-14 08:53:15 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4KaonMinusAbsorption physics process ------
//                   by Larry Felawka (TRIUMF), April 1998
//                     E-mail: felawka@alph04.triumf.ca
// ************************************************************
//-----------------------------------------------------------------------------

#ifndef G4KaonMinusAbsorption_h
#define G4KaonMinusAbsorption_h 1
 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4GHEKinematicsVector.hh"
// Class Description
// Alternative process for absorption of K- at rest; 
// to be used in your physics list in case you need this physics.
// Class Description - End


class G4KaonMinusAbsorption : public G4VRestProcess
 
{ 
  private:
  // hide assignment operator as private 
      G4KaonMinusAbsorption& operator=(const G4KaonMinusAbsorption &right);
      G4KaonMinusAbsorption(const G4KaonMinusAbsorption& );
   
  public:
 
     G4KaonMinusAbsorption(const G4String& processName ="KaonMinusAbsorption");
 
    ~G4KaonMinusAbsorption();

     G4bool IsApplicable(const G4ParticleDefinition&);

  // null physics table
     void BuildPhysicsTable(const G4ParticleDefinition&){}

     G4double AtRestGetPhysicalInteractionLength(const G4Track&,
						 G4ForceCondition*);

  // zero mean lifetime
     G4double GetMeanLifeTime(const G4Track& aTrack,
			      G4ForceCondition* condition) {return 0.0;}

     G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

  // return number of secondaries produced
     G4int GetNumberOfSecondaries();

  // pointer to array containg kinematics of secondaries
     G4GHEKinematicsVector* GetSecondaryKinematics();

  private:

     void GenerateSecondaries();
     void Poisso( G4float, G4int* );
     void Normal( G4float* );
     void KaonMinusAbsorption( G4int* );
     G4int NFac( G4int );

  private:

// global time-of-flight of stopped AntiProton
     G4float  globalTime;

// atomic mass of target nucleus
     G4float  targetAtomicMass;

// charge of target nucleus
     G4float  targetCharge;

     G4GHEKinematicsVector* pv;
     G4GHEKinematicsVector* eve;
     G4GHEKinematicsVector* gkin;

     G4float  evapEnergy1;
     G4float  evapEnergy3;

     G4int    ngkine;

     G4int    ntot;
     G4GHEKinematicsVector result;

     G4float  massKaonMinus;
     G4float  massGamma;
     G4float  massPionZero;
     G4float  massProton;
     G4float  massLambda;

     G4ParticleDefinition* pdefKaonMinus;
     G4ParticleDefinition* pdefGamma;
     G4ParticleDefinition* pdefPionZero;
     G4ParticleDefinition* pdefProton;
     G4ParticleDefinition* pdefNeutron;
     G4ParticleDefinition* pdefLambda;
     G4ParticleDefinition* pdefDeuteron;
     G4ParticleDefinition* pdefTriton;
     G4ParticleDefinition* pdefAlpha;

};

#endif
 
