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
// $Id: G4MuMinusCapturePrecompound.hh 79691 2014-03-12 12:58:09Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class header file 
//
// File name:  G4MuMinusCapturePrecompound
//
// Author:        V.Ivanchenko 
// 
// Creation date: 24 April 2012 on base of G4MuonMinusCaptureAtRest
//
// Class Description: 
//
// Sampling of mu- capture by atomic nucleus
//
//-----------------------------------------------------------------------------
//
// Modifications: 
//
//-----------------------------------------------------------------------------

#ifndef G4MuMinusCapturePrecompound_h
#define G4MuMinusCapturePrecompound_h 1

#include "globals.hh"
#include "G4Nucleus.hh"
#include "G4Track.hh"
#include "G4HadProjectile.hh"
#include "G4HadSecondary.hh"
#include "G4HadFinalState.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4Fancy3DNucleus.hh"

class G4VPreCompoundModel;

class G4MuMinusCapturePrecompound : public G4HadronicInteraction
{ 
public:
 
  G4MuMinusCapturePrecompound(G4VPreCompoundModel* ptr=0);
 
  ~G4MuMinusCapturePrecompound();

  G4HadFinalState* ApplyYourself(const G4HadProjectile &aTrack, 
				 G4Nucleus & targetNucleus );

  void ModelDescription(std::ostream& outFile) const; 

private:

  inline void AddNewParticle(const G4ParticleDefinition* aParticle,
			     G4ThreeVector& direction,
			     G4double kinEnergy);

  // hide assignment operator as private 
  G4MuMinusCapturePrecompound& operator=(const G4MuMinusCapturePrecompound &right);
  G4MuMinusCapturePrecompound(const G4MuMinusCapturePrecompound& );

  G4HadFinalState result;
  G4Fancy3DNucleus fNucleus;
  const G4ParticleDefinition* fProton;
  const G4ParticleDefinition* fNeutron;
  G4VPreCompoundModel* fPreCompound;
  G4double fMuMass;
  G4double fThreshold;
  G4double fTime;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4MuMinusCapturePrecompound::AddNewParticle(
                      const G4ParticleDefinition* aParticle,
		      G4ThreeVector& direction,
		      G4double kinEnergy)
{
  G4DynamicParticle* dp = new G4DynamicParticle(aParticle,
                                                direction,
                                                kinEnergy);
  G4HadSecondary hs(dp);
  hs.SetTime(fTime);
  result.AddSecondary(hs);
}

#endif


