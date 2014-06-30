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
// G4 Low energy model: n-n or p-p scattering
// F.W. Jones, L.G. Greeniaus, H.P. Wellisch
//  
// For further comments see G4LEppData.hh and G4LEpp.cc
//
// 30.01.14 V. Grichine add SampleInvariantT and inherit 
//             from G4HadronElastic

#ifndef G4LEpp_h
#define G4LEpp_h 1
 
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronElastic.hh"

class G4LEpp : public G4HadronElastic
{
private:

  enum { NENERGY=40, NANGLE=180 };

public:

  G4LEpp();

  virtual ~G4LEpp();
 
  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
  				 G4Nucleus& targetNucleus);

  G4double SampleInvariantT(const G4ParticleDefinition* p, 
			    G4double plab, G4int Z, G4int A);
  
private:

  // The following arrays are declared static to allow the use of initializers.
  // They are initialized in G4LEppData.hh

  // Coulomb effects suppressed:
  static const G4float Sig[NENERGY][NANGLE];
  static const G4float elab[NENERGY]; 
  static const G4float dSigmax[NENERGY];
  static const G4float Sigtot[NENERGY];

};

#endif
