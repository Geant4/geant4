//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: HadrontherapyPhysicsList.hh,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#ifndef HadrontherapyPhysicsList_h
#define HadrontherapyPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class HadrontherapyPhysicsListMessenger;


class HadrontherapyPhysicsList: public G4VModularPhysicsList
{
public:
  HadrontherapyPhysicsList();
  virtual ~HadrontherapyPhysicsList();

  virtual void SetCuts();
  void AddPhysicsList(const G4String& name);  

  // Production thresholds, expressed in range
  void SetGammaCut(G4double cut);
  void SetElectronCut(G4double cut);
  void SetProtonCut(G4double cut);
  
private:
  G4bool electronIsRegistered;
  G4bool positronIsRegistered;
  G4bool photonIsRegistered;
  G4bool ionIsRegistered;
  G4bool protonHadroIsRegistered;
  G4bool chargedParticleIsRegistered;
  G4bool muonIsRegistered;
  G4bool decayIsRegistered;

  G4double currentDefaultCut;
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;

  HadrontherapyPhysicsListMessenger* messenger;
};

#endif



