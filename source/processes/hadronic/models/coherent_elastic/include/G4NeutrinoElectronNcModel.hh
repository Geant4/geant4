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
// $Id: G4NeutrinoElectronNcModel.hh 90228 2015-05-21 08:49:57Z gcosmo $
//
// Geant4 Header : G4NeutrinoElectronNcModel
//
// Author : V.Grichine 6.4.17
//  
// Modified:
//
// Class Description
// Default model for neutrino-electron elastic (neutral current) scattering; 
// Class Description - End

#ifndef G4NeutrinoElectronNcModel_h
#define G4NeutrinoElectronNcModel_h 1
 
#include "globals.hh"
#include "G4HadronElastic.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"

class G4ParticleDefinition;

class G4NeutrinoElectronNcModel : public G4HadronElastic
{
public:

  G4NeutrinoElectronNcModel(const G4String& name = "nu-e-elastic");

  virtual ~G4NeutrinoElectronNcModel();

  virtual G4bool IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus);

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
					  G4Nucleus & targetNucleus);

  // sample recoil electron energy in lab frame

  G4double SampleElectronTkin(const G4HadProjectile* aParticle);

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};


  
  virtual void ModelDescription(std::ostream&) const;

private:

  G4ParticleDefinition* theElectron; 
  G4double fSin2tW;    // sin^2theta_Weinberg
  G4double fCutEnergy; // minimal recoil electron energy detected

};



#endif
