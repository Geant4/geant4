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
// G4 Low energy model: n-p scattering
// F.W. Jones, L.G. Greeniaus, H.P. Wellisch
//  
// For further comments see G4LEnpData.hh and G4LEnp.cc
//
// 02.01.2014 V. Grichine add SampInvariantT for testing ds/dt

#ifndef G4LEnp_h
#define G4LEnp_h 1
 
#include "globals.hh"
#include "G4HadronElastic.hh"

class G4LEnp : public G4HadronElastic  // G4HadronicInteraction
{
private:

  enum { NENERGY=39, NANGLE=180 };

public:

  explicit G4LEnp();

  ~G4LEnp() override;
 
  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                 G4Nucleus& targetNucleus) override;

  // sample momentum transfer using Lab. momentum

  G4double SampleInvariantT(const G4ParticleDefinition* p, 
			    G4double plab, G4int Z, G4int A) override;
private:

 // The following arrays are declared static to allow the use of initializers.
 // They are initialized in G4LEnpData.hh

  static const G4float sig[NENERGY][NANGLE];
  static const G4float pcm[NENERGY], elab[NENERGY]; 
  static const G4float dsigmax[NENERGY], sigtot[NENERGY];

};

#endif
