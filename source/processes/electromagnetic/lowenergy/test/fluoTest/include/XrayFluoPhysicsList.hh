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
#ifndef XrayFluoPhysicsList_h
#define XrayFluoPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class XrayFluoPhysicsListMessenger;
class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;
class XrayFluoDetectorConstruction;
class XrayFluoPhysicsList: public G4VUserPhysicsList
{
public:

  XrayFluoPhysicsList(XrayFluoDetectorConstruction*);
  virtual ~XrayFluoPhysicsList();
  
protected:

  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
protected:

  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  virtual void ConstructBarions();
  virtual void ConstructIons();

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
  void SetProtonCut(G4double);
  void SetCutsByEnergy(G4double);
  

private:

  G4LowEnergyIonisation*  LeIoprocess;
  G4LowEnergyPhotoElectric* LePeprocess;
  G4LowEnergyBremsstrahlung* LeBrprocess;
  XrayFluoPhysicsListMessenger* physicsListMessenger;
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;
  XrayFluoDetectorConstruction* pDet;
  
};

#endif







