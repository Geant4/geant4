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
// $Id: ExGflashDetectorConstruction.hh 98248 2016-07-04 17:25:20Z gcosmo $
//
/// \file parameterisations/gflash/include/ExGflashDetectorConstruction.hh
/// \brief Definition of the ExGflashDetectorConstruction class
//
#ifndef ExGflashDetectorConstruction_h
#define ExGflashDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Region;

class GFlashHomoShowerParameterisation;
class GFlashShowerModel;
class GFlashHitMaker;
class GFlashParticleBounds;

class ExGflashDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  ExGflashDetectorConstruction();
  ~ExGflashDetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  
  const G4VPhysicalVolume* GetCristal(int num__crystal)
  {return fCrystal_phys[num__crystal];};
  
  
private:
  G4LogicalVolume*    fCrystal_log;
  G4VPhysicalVolume*  fCrystal_phys[100];
  G4Region*           fRegion;

  static G4ThreadLocal GFlashShowerModel* fFastShowerModel; 
  static G4ThreadLocal GFlashHomoShowerParameterisation* fParameterisation;
  static G4ThreadLocal GFlashParticleBounds* fParticleBounds;
  static G4ThreadLocal GFlashHitMaker* fHitMaker;
};

#endif




















