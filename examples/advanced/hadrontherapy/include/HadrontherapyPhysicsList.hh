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

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class HadrontherapyDetectorConstruction;
class HadrontherapyPhysicsListMessenger;

// ----------------------------------------------------------------
class HadrontherapyPhysicsList: public G4VUserPhysicsList
{
public:
  HadrontherapyPhysicsList( HadrontherapyDetectorConstruction*);
  ~HadrontherapyPhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBarions();
  void ConstructIons();

protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructHad();
  void ConstructOp();

public:  
  void SetGammaCut(G4double);
  void SetElectronCut(G4double);
  void SetProtonCut(G4double);
  void GetRange(G4double);
  void SetMaxStep(G4double);
  
private:
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;
  G4double currentDefaultCut;
  G4double defaultCutValue;
  HadrontherapyDetectorConstruction* pDet;
  HadrontherapyPhysicsListMessenger* physicsListMessenger;
};

#endif



