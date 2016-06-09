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
// $Id: XrayFluoPhysicsList.hh
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoPhysicsList_h
#define XrayFluoPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
//#include "G4LivermorePhotoElectricModel.hh"

//#include "XrayFluoPlaneDetectorConstruction.hh"
//#include "XrayFluoDetectorConstruction.hh"
//#include "XrayFluoMercuryConstruction.hh"


class G4LowEnergyIonisation;
class G4LivermorePhotoElectricModel;
class G4LowEnergyBremsstrahlung;
class G4eIonisation;


class XrayFluoPhysicsListMessenger;
class XrayFluoDetectorConstruction;
class XrayFluoPlaneDetectorConstruction;
class XrayFluoMercuryDetectorConstruction;


class XrayFluoPhysicsList: public G4VUserPhysicsList
{
public:

  XrayFluoPhysicsList(XrayFluoDetectorConstruction*);
  XrayFluoPhysicsList(XrayFluoPlaneDetectorConstruction*);
  XrayFluoPhysicsList(XrayFluoMercuryDetectorConstruction*);


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

public:
 
  void SetCuts();
  void SetGammaCut(G4double);
  void SetElectronCut(G4double);

  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);

  void SetGELowLimit(G4double);

  void SetLowEnSecPhotCut(G4double);
//  void SetLowEnSecElecCut(G4double);
  void SetProtonCut(G4double);
  void SetCutsByEnergy(G4double);
  

private:


  G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel;	
  XrayFluoPhysicsListMessenger* physicsListMessenger;
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;
  XrayFluoDetectorConstruction* pDet;
  XrayFluoPlaneDetectorConstruction* planeDet;
  XrayFluoMercuryDetectorConstruction* mercuryDet;

};

#endif







