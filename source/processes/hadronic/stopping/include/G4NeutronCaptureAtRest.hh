// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronCaptureAtRest.hh,v 1.3 2000-12-14 08:53:15 hpw Exp $
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
//      ------------ G4NeutronCaptureAtRest physics process ------
//                   by Larry Felawka (TRIUMF), April 1998
//                     E-mail: felawka@alph04.triumf.ca
// ************************************************************
//-----------------------------------------------------------------------------

#ifndef G4NeutronCaptureAtRest_h
#define G4NeutronCaptureAtRest_h 1
// Class Description
// Process for capture of neutrons at rest; 
// to be used in your physics list in case you need this physics.
// Class Description - End

 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4GHEKinematicsVector.hh"

class G4NeutronCaptureAtRest : public G4VRestProcess
 
{ 
  private:
  // hide assignment operator as private 
      G4NeutronCaptureAtRest& operator=(const G4NeutronCaptureAtRest &right);
      G4NeutronCaptureAtRest(const G4NeutronCaptureAtRest& );
   
  public:
 
     G4NeutronCaptureAtRest(const G4String& processName ="NeutronCaptureAtRest");
 
    ~G4NeutronCaptureAtRest();

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
     void Normal( G4float* );
     void NeutronCapture( G4int* );
     G4double AtomAs( G4float, G4float );

  private:

// global time-of-flight of stopped hadron
     G4float  globalTime;

// atomic mass of target nucleus
     G4float  targetAtomicMass;

// charge of target nucleus
     G4float  targetCharge;

     G4GHEKinematicsVector* pv;
     G4GHEKinematicsVector* eve;
     G4GHEKinematicsVector* gkin;

     G4int    ngkine;

     G4int    ntot;
     G4GHEKinematicsVector result;

     G4float  massProton;
     G4float  massNeutron;
     G4float  massElectron;
     G4float  massDeuteron;
     G4float  massAlpha;

     G4ParticleDefinition* pdefGamma;
     G4ParticleDefinition* pdefNeutron;

};

#endif
 
