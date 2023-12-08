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
// Geant4 header : G4NeutronRadCaptureHP
// Created:  02 October 2023
// Author  V.Ivanchenko
//  
// Modified:
//
// Class Description
// Sampling of neutron radiative capture. The same approach is used
// as in the G4NeutronRadCapture with the difference in sampling of
// nucleus termal motion.
// Class Description - End
//

#ifndef G4NeutronRadCaptureHP_h
#define G4NeutronRadCaptureHP_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

class G4VEvaporationChannel;
class G4IonTable;
class G4ParticleHPManager;

class G4NeutronRadCaptureHP : public G4HadronicInteraction
{
public:

  G4NeutronRadCaptureHP();

  ~G4NeutronRadCaptureHP() override;

  G4HadFinalState* ApplyYourself(const G4HadProjectile & aTrack, 
                                 G4Nucleus & targetNucleus) override;

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  G4NeutronRadCaptureHP & operator=
  (const G4NeutronRadCaptureHP &right) = delete;
  G4NeutronRadCaptureHP(const G4NeutronRadCaptureHP&) = delete;

private:

  G4int icID{-1}; // creator model ID for e- produced by internal conversion
  G4int secID{-1};  // creator model ID for the other secondaries produced by this model
  const G4ParticleDefinition* electron;
  G4ParticleHPManager* fManagerHP;
  G4VEvaporationChannel* photonEvaporation{nullptr};
  G4IonTable* theTableOfIons;
  G4double lowestEnergyLimit;
  G4double minExcitation;
  G4double emax;
  G4double emaxT;
  G4LorentzVector lab4mom;
  G4bool fLocalPE{false};

};

#endif
