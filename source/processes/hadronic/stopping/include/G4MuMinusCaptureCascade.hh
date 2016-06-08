// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuMinusCaptureCascade physics process --------
//                   by Vladimir Ivanchenko
//                     E-mail: Vladimir.Ivantchenko@cern.ch
//                            April 2000
// **************************************************************
//-----------------------------------------------------------------------------

#ifndef G4MuMinusCaptureCascade_h
#define G4MuMinusCaptureCascade_h 1
#include "g4std/iomanip" 
#include "globals.hh"
#include "Randomize.hh" 

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4GHEKinematicsVector.hh"
#include "G4ParticleTypes.hh"
#include "G4DynamicParticle.hh"

class G4MuMinusCaptureCascade
 
{ 
  private:
  // hide assignment operator as private 
      G4MuMinusCaptureCascade& operator=(const G4MuMinusCaptureCascade &right);
      G4MuMinusCaptureCascade(const G4MuMinusCaptureCascade& );
   
  public:
 
      G4MuMinusCaptureCascade();
 
     ~G4MuMinusCaptureCascade();

      G4int DoCascade(const G4double Z, const G4double massA, 
                            G4GHEKinematicsVector* Cascade);

      void DoBoundMuonMinusDecay(G4double Z, G4double massA, 
                                 G4int* nCascade, G4GHEKinematicsVector* Cascade);

  private:

      G4double GetKShellEnergy(G4double Z);

      G4double GetLinApprox(const size_t N, const G4double X[], const G4double Y[], 
                            G4double Xuser);

      G4ThreeVector GetRandomVec();


      void AddNewParticle(G4ParticleDefinition* aParticle,
                          G4ThreeVector Momentum,
                          G4double mass,
                          G4int* nParticle,
                          G4GHEKinematicsVector* Cascade);

      G4double Emass, MuMass;
      G4ParticleDefinition* theElectron;
      G4ParticleDefinition* theGamma;

};

#endif


