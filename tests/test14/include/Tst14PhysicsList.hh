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
// $Id: Tst14PhysicsList.hh,v 1.15 2010-08-29 19:50:17 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Original author unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Redesigned for modular PhysicsList
// 06 Nov 2003 MGP          Introduced test of Bremsstrahlung angular distributions
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

private:

  G4bool electronIsRegistered;
  G4bool positronIsRegistered;
  G4bool photonIsRegistered;

  G4double cutForGamma;
  G4double cutForElectron;

  Tst14PhysicsListMessenger* messenger;

};

#endif







