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
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuMinusCaptureCascade physics process --------
//                   by Vladimir Ivanchenko
//                     E-mail: Vladimir.Ivantchenko@cern.ch
//                            April 2000
//
//-----------------------------------------------------------------------------

#ifndef G4MuMinusCaptureCascade_h
#define G4MuMinusCaptureCascade_h 1

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

private:

  // hide assignment operator as private 
  G4MuMinusCaptureCascade& operator=(const G4MuMinusCaptureCascade &right);
  G4MuMinusCaptureCascade(const G4MuMinusCaptureCascade& );

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


