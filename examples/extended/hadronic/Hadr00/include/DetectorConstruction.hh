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
/// \file hadronic/Hadr00/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 109182 2018-04-03 07:18:45Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorConstruction
//
// Created: 20.06.2008 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
// 

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class DetectorMessenger;
class G4VModularPhysicsList;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  void SetWorldMaterial(const G4String&);
  void SetTargetMaterial(const G4String&);

  void SetTargetRadius(G4double val);
  void SetTargetLength(G4double val);

  const G4Material* GetTargetMaterial() const { return fTargetMaterial; }
  G4VModularPhysicsList* GetPhysicsList() { return fPhysList; }
  void SetPhysicsList(G4VModularPhysicsList* ptr) { fPhysList = ptr; }

private:

  void ComputeGeomParameters();

  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);

  G4double fRadius;
  G4double fLength;
  G4double fWorldR;
  G4double fWorldZ;

  G4Material*  fTargetMaterial;
  G4Material*  fWorldMaterial;

  G4Tubs* fSolidW;
  G4Tubs* fSolidA;

  G4LogicalVolume* fLogicTarget;
  G4LogicalVolume* fLogicWorld;

  G4VPhysicalVolume* fPhysWorld;

  DetectorMessenger* fDetectorMessenger;
  G4VModularPhysicsList* fPhysList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif

