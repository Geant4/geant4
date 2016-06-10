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
// $Id: ExGflashDetectorConstruction.hh 73006 2013-08-15 08:17:11Z gcosmo $
//
/// \file parameterisations/gflash/include/ExGflashDetectorConstruction.hh
/// \brief Definition of the ExGflashDetectorConstruction class
//
#ifndef ExGflashDetectorConstruction_h
#define ExGflashDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class GFlashHomoShowerParameterisation;
class GFlashHitMaker;
class GFlashShowerModel;
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
  G4LogicalVolume* fExperimentalHall_log;
  G4LogicalVolume* fCalo_log;
  
  G4VPhysicalVolume* fExperimentalHall_phys;      
  G4VPhysicalVolume* fCalo_phys;
  
  G4Box *fExperimentalHall_box;
  
  G4double fExperimentalHall_x;
  G4double fExperimentalHall_y;
  G4double fExperimentalHall_z;
  
  G4double fCalo_xside;
  G4double fCalo_yside;
  G4double fCalo_zside;    
  
  G4int    fNbOfCrystals;                // Nb of chambers in the tracker region
  G4double fCrystalWidth;                // width of the chambers
  G4double fCrystalLenght;
  //@@@  ExGflashDetectorConstruction : wie mache ich das am besten ?
  G4Box *fCrystal[100];      
  G4LogicalVolume* fCrystal_log[100];
  G4VPhysicalVolume*  fCrystal_phys[100];
  
  // Gflash members    
  GFlashHomoShowerParameterisation* fTheParameterisation;
  GFlashHitMaker*                   fTheHMaker;
  GFlashParticleBounds*             fTheParticleBounds;
  GFlashShowerModel*                fTheFastShowerModel;     
  G4Region*                         fRegion;
};

#endif




















