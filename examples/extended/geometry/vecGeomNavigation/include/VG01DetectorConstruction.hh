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
/// \file VG01DetectorConstruction.hh
/// \brief Definition of the VG01DetectorConstruction class


#ifndef VG01DetectorConstruction_h
#define VG01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4GDMLParser.hh"
#include "G4String.hh"

class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class VG01DetectorMessenger;

class VG01DetectorConstruction : public G4VUserDetectorConstruction {

public:

  VG01DetectorConstruction();
  ~VG01DetectorConstruction();

  G4VPhysicalVolume* Construct() override;

  void SetGDMLFileName ( const G4String& gdmlfile ) { fGDMLFileName = gdmlfile; }
  void SetUseVecGeom (bool b) { fUseVecGeom = b; }
  void SetMagFieldValue(const G4double fieldValue ) { fglobFieldValue = fieldValue; }

  static G4double GetFieldValue() { return fglobFieldValue; }

private:
  void CreateMagFieldAndIntegrator();

private:
  // this static member is for the print out
  static G4double        fglobFieldValue;

  G4String               fGDMLFileName;
  G4GDMLParser           fParser;
  G4VPhysicalVolume*     fWorld;
  G4FieldManager*        fFieldMgr;
  G4UniformMagField*     fUniformMagField;
  bool                   fUseVecGeom = true;
  VG01DetectorMessenger* fDetectorMessenger= nullptr;
};

#endif
