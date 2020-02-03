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
// Geant4 Header : G4NeutrinoElectronCcModel
//
// Author : V.Grichine 26.4.17
//  
// Modified:
//
// Class Description
// Default model for neutrino-electron 'inelastic' (charge current) scattering; 
// Class Description - End

#ifndef G4NeutrinoElectronCcModel_h
#define G4NeutrinoElectronCcModel_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"

class G4ParticleDefinition;

class G4NeutrinoElectronCcModel : public G4HadronicInteraction
{
public:

  G4NeutrinoElectronCcModel(const G4String& name = "nu-e-elastic");

  virtual ~G4NeutrinoElectronCcModel();

  virtual G4bool IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus);

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
					  G4Nucleus & targetNucleus);

  // sample recoil electron energy in lab frame

  G4double SampleCosCMS(const G4HadProjectile* aParticle);

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};


  
  virtual void ModelDescription(std::ostream&) const;

private:

  G4ParticleDefinition* theNeutrinoE;
  G4ParticleDefinition* theAntiNeutrinoE;
  G4ParticleDefinition* theNeutrinoMu;
  G4ParticleDefinition* theAntiNeutrinoMu;
  G4ParticleDefinition* theNeutrinoTau;
  G4ParticleDefinition* theAntiNeutrinoTau;
  
  G4ParticleDefinition* theMuonMinus;
  G4ParticleDefinition* theTauMinus;
 
  G4double fSin2tW;    // sin^2theta_Weinberg
  G4double fCutEnergy; // minimal recoil electron energy detected

};



#endif
