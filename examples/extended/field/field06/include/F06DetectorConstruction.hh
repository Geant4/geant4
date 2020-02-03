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
/// \file field/field06/include/F06DetectorConstruction.hh
/// \brief Definition of the F06DetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F06DetectorConstruction_h
#define F06DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
 
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

class G4UniformGravityField;
class G4FieldManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F06DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
 
    F06DetectorConstruction();
    virtual ~F06DetectorConstruction();

  public:
 
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:
 
     G4Material*        fVacuum;

     G4Box*             fSolidWorld;    //pointer to the solid World
     G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
     G4VPhysicalVolume* fPhysiWorld;    //pointer to the physical World

     static G4ThreadLocal G4UniformGravityField* fField;

  private:
 
     void DefineMaterials();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
