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
// $Id: G4EmCaptureCascade.hh 101422 2016-11-17 10:41:23Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class header file 
//
// File name:  G4EmCaptureCascade
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 22 April 2012 on base of G4MuMinusCaptureCascade
//
// Class Description: 
//
// Simulation of electromagnetic cascade from capture level to K-shell
// of the mesonic atom
//
// Probabilities of gamma and Auger transitions from
// N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.
//
//-----------------------------------------------------------------------------
//
// Modifications: 
//
//-----------------------------------------------------------------------------

#ifndef G4EmCaptureCascade_h
#define G4EmCaptureCascade_h 1

#include "globals.hh"
#include "G4Nucleus.hh"
#include "G4Track.hh"
#include "G4HadProjectile.hh"
#include "G4HadSecondary.hh"
#include "G4HadFinalState.hh"
#include "G4HadronicInteraction.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"

class G4EmCaptureCascade : public G4HadronicInteraction
{ 
public:
 
  explicit G4EmCaptureCascade();
 
  virtual ~G4EmCaptureCascade();

  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile &aTrack, 
					 G4Nucleus & targetNucleus );

  virtual void ModelDescription(std::ostream& outFile) const; 

private:

  inline void AddNewParticle(G4ParticleDefinition* aParticle,
			     G4double kinEnergy);

  // hide assignment operator as private 
  G4EmCaptureCascade& operator=(const G4EmCaptureCascade &right) = delete;
  G4EmCaptureCascade(const G4EmCaptureCascade& ) = delete;

  G4HadFinalState result;
  G4ParticleDefinition* theElectron;
  G4ParticleDefinition* theGamma;
  G4double fMuMass;
  G4double fTime;
  G4double fLevelEnergy[14];
  G4double fKLevelEnergy[93];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4EmCaptureCascade::AddNewParticle(G4ParticleDefinition* aParticle,
				   G4double kinEnergy)
{
  G4DynamicParticle* dp = new G4DynamicParticle(aParticle,
                                                G4RandomDirection(),
                                                kinEnergy);
  G4HadSecondary hs(dp);
  hs.SetTime(fTime);
  result.AddSecondary(hs);
}

#endif


