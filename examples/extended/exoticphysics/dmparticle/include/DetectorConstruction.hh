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
/// \file exoticphysics/dmparticle/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 68036 2013-03-13 14:13:45Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  // set geometry parameters
  void SetSizeZ(G4double);
  void SetSizeXY(G4double);              
  void SetMaterial(const G4String&);            
          
  // access to geometry
  inline G4double     GetWorldSizeZ()    {return fWorldSizeZ;};
  inline G4double     GetAbsorSizeZ()    {return fAbsorSizeZ;};
  inline const G4Material* GetAbsorMaterial() {return fAbsorMaterial;};
                           
private:

  void PrintParameters();
  
  G4double            fWorldSizeZ;
  G4double            fWorldSizeXY;
  G4Material*         fWorldMaterial;           
  G4double            fAbsorSizeZ;
  G4double            fAbsorSizeXY;
  G4Material*         fAbsorMaterial;

  G4LogicalVolume*    fLogAbsor;
  G4VPhysicalVolume*  fWorld;
               
  DetectorMessenger*  fDetectorMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

