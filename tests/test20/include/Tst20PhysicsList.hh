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
#ifndef Tst20PhysicsList_h
#define Tst20PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst20PhysicsListMessenger;
class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class Tst20PhysicsList: public G4VUserPhysicsList
{
public:

  Tst20PhysicsList();
  virtual ~Tst20PhysicsList();
  
protected:

  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
protected:

  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  
protected:
  // these methods Construct physics processes and register them
  void ConstructEM();
  void ConstructGeneral();

  void SetCuts();

public:

  void SetGammaCut(G4double);
  void SetElectronCut(G4double);

  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);
  void SetGELowLimit(G4double);
  void SetLowEnSecPhotCut(G4double);
  void SetLowEnSecElecCut(G4double);

private:

  G4LowEnergyIonisation*  lowEIoniProcess;
  G4LowEnergyPhotoElectric* lowEPhotoelProcess;
  G4LowEnergyBremsstrahlung* lowEBremProcess;
  Tst20PhysicsListMessenger* physicsListMessenger;
  G4double cutForGamma;
  G4double cutForElectron;

};

#endif







