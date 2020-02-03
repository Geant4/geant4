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
/// \file ExGflash2DetectorConstruction.hh
/// \brief Definition of the ExGflash2DetectorConstruction class
//
#ifndef ExGflash2DetectorConstruction_h
#define ExGflash2DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class GFlashHomoShowerParameterisation;
class GFlashShowerModel;
class GFlashHitMaker;
class GFlashParticleBounds;

class ExGflash2DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  ExGflash2DetectorConstruction();
  ~ExGflash2DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  
  const G4VPhysicalVolume* GetCristal(int aNumCrystal)
  {return fCrystalPhys[aNumCrystal];};
  
  
private:
  G4LogicalVolume*    fCrystalLog;
  G4VPhysicalVolume*  fCrystalPhys[100];

  static G4ThreadLocal GFlashShowerModel* fFastShowerModel; 
  static G4ThreadLocal GFlashHomoShowerParameterisation* fParameterisation;
  static G4ThreadLocal GFlashParticleBounds* fParticleBounds;
  static G4ThreadLocal GFlashHitMaker* fHitMaker;
};

#endif




















