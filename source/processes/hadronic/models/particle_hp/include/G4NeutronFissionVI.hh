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
// Geant4 header : G4NeutronFissionVI
// Created:  03 October 2023
// Author  V.Ivanchenko
//  
// Modified:
//
// Class Description
// Sampling of neutron induced fission using de-excitation module
// Class Description - End
//

#ifndef G4NeutronFissionVI_h
#define G4NeutronFissionVI_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

class G4VEvaporationChannel;
class G4IonTable;
class G4ExcitationHandler;
class G4ParticleHPManager;

class G4NeutronFissionVI : public G4HadronicInteraction
{
public:

  G4NeutronFissionVI();

  ~G4NeutronFissionVI() override;

  G4HadFinalState* ApplyYourself(const G4HadProjectile & aTrack, 
                                 G4Nucleus & targetNucleus) override;

  void InitialiseModel() override;

  G4NeutronFissionVI & operator=(const G4NeutronFissionVI &right) = delete;
  G4NeutronFissionVI(const G4NeutronFissionVI&) = delete;

private:

  G4int secID{-1};  // creator model ID for secondaries produced by this model
  G4ParticleHPManager* fManagerHP;
  G4ExcitationHandler* fHandler{nullptr};
  G4VEvaporationChannel* fFission{nullptr};
  G4IonTable* theTableOfIons;

  G4double minExcitation;
  G4double emaxT;
  G4LorentzVector lab4mom;
  G4bool fLocalHandler{false};
};

#endif
