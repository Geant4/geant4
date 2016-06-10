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
//     $Id: G4AntiNeutronAnnihilationAtRest.hh 66367 2012-12-18 09:18:08Z gcosmo $
//     G4AntiNeutronAnnihilationAtRest physics process
//     Larry Felawka (TRIUMF), April 1998
//---------------------------------------------------------------------

#ifndef G4AntiNeutronAnnihilationAtRest_h
#define G4AntiNeutronAnnihilationAtRest_h 1

// Class Description:
//
// Process for absorption of n-bar at rest. 
// To be used in your physics list in case you need this physics.
 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4GHEKinematicsVector.hh"
#include "G4HadronicProcessType.hh"

class G4AntiNeutronAnnihilationAtRest : public G4VRestProcess
 
{ 
  private:
  // hide assignment operator as private 
      G4AntiNeutronAnnihilationAtRest& operator=(const G4AntiNeutronAnnihilationAtRest &right);
      G4AntiNeutronAnnihilationAtRest(const G4AntiNeutronAnnihilationAtRest& );
   
  public:
 
     G4AntiNeutronAnnihilationAtRest(const G4String& processName ="AntiNeutronAnnihilationAtRest", 
                       G4ProcessType   aType = fHadronic );
 
    ~G4AntiNeutronAnnihilationAtRest();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void PreparePhysicsTable(const G4ParticleDefinition&);

     void BuildPhysicsTable(const G4ParticleDefinition&);

     G4double AtRestGetPhysicalInteractionLength(const G4Track&,
						 G4ForceCondition*);

  // zero mean lifetime
     G4double GetMeanLifeTime(const G4Track& ,
			      G4ForceCondition* ) {return 0.0;}

     G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

  // return number of secondaries produced
     G4int GetNumberOfSecondaries();

  // pointer to array containg kinematics of secondaries
     G4GHEKinematicsVector* GetSecondaryKinematics();

  private:

     void GenerateSecondaries();
     void Poisso( G4float, G4int* );
     void Normal( G4float* );
     void AntiNeutronAnnihilation( G4int* );
     G4double ExNu( G4float );
     G4int NFac( G4int );

  private:

// global time-of-flight of stopped AntiNeutron
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

     G4float  massPionMinus;
     G4float  massPionZero;
     G4float  massPionPlus;
     G4float  massGamma;
     G4float  massAntiNeutron;
     G4float  massNeutron;

     G4ParticleDefinition* pdefGamma;
     G4ParticleDefinition* pdefPionPlus;
     G4ParticleDefinition* pdefPionZero;
     G4ParticleDefinition* pdefPionMinus;
     G4ParticleDefinition* pdefProton;
     G4ParticleDefinition* pdefNeutron;
     G4ParticleDefinition* pdefAntiNeutron;
     G4ParticleDefinition* pdefDeuteron;
     G4ParticleDefinition* pdefTriton;
     G4ParticleDefinition* pdefAlpha;

};

#endif
 
