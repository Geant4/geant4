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
// $Id: Tst14PhysicsList.hh,v 1.12 2003-02-23 17:22:07 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Original author unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Redesigned for modular PhysicsList
//
// -------------------------------------------------------------------

// Class description:
// System test for e/gamma, standard photon processes for PhysicsList
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef TST14PHYSICSLIST_HH
#define TST14PHYSICSLIST_HH 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class Tst14PhysicsListMessenger;

class Tst14PhysicsList: public G4VModularPhysicsList {
public:
  
  Tst14PhysicsList();

  virtual ~Tst14PhysicsList();

  virtual void SetCuts();
  
  // Register PhysicsList chunks
  void AddPhysicsList(const G4String& name);

  // Production thresholds, expressed in range
  void SetGammaCut(G4double cut);
  void SetElectronCut(G4double cut);

  // Production thresholds, expressed in energy, for photons, electrons and both
  void SetGammaLowLimit(G4double cut);
  void SetElectronLowLimit(G4double cut);
  void SetGELowLimit(G4double cut);

  // Cut for generation of secondaries for EEDL/EPDL processes
  void SetLowEnSecPhotCut(G4double cut);
  void SetLowEnSecElecCut(G4double cut);
 
  // Activation of Auger effect in electron ionisation and photoelectric effect
  void ActivateAuger(G4bool flag);

private:

  G4bool electronIsRegistered;
  G4bool positronIsRegistered;
  G4bool photonIsRegistered;

  G4double cutForGamma;
  G4double cutForElectron;

  Tst14PhysicsListMessenger* messenger;

};

#endif







