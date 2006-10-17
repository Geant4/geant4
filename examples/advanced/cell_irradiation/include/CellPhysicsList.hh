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
//    **************************************
//    *                                    *
//    *         CellPhysicsList.hh         *
//    *                                    *
//    **************************************
//
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

#ifndef CELLPHYSICSLIST_HH
#define CELLPHYSICSLIST_HH 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class CellPhysicsListMessenger;

class CellPhysicsList: public G4VModularPhysicsList {
public:
  
  CellPhysicsList();

  virtual ~CellPhysicsList();

  virtual void SetCuts();

  void AddPhysicsList(const G4String& name);
  // Register PhysicsList chunks

  void SetParticleCut(G4double value);
  // Production thresholds, expressed in range

private:

  G4bool electronIsRegistered;
  G4bool positronIsRegistered;
  G4bool photonIsRegistered;
  G4bool protonIsRegistered;
  
  CellPhysicsListMessenger* messenger;

};
#endif







