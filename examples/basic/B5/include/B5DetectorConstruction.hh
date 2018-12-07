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
/// \file B5DetectorConstruction.hh
/// \brief Definition of the B5DetectorConstruction class

#ifndef B5DetectorConstruction_h
#define B5DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"

#include <vector>

class B5MagneticField;

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

/// Detector construction

class B5DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B5DetectorConstruction();
    virtual ~B5DetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void SetArmAngle(G4double val);
    G4double GetArmAngle() { return fArmAngle; }
    
    void ConstructMaterials();
    
  private:
    void DefineCommands();

    G4GenericMessenger* fMessenger;
    
    static G4ThreadLocal B5MagneticField* fMagneticField;
    static G4ThreadLocal G4FieldManager* fFieldMgr;
    
    G4LogicalVolume* fHodoscope1Logical;
    G4LogicalVolume* fHodoscope2Logical;
    G4LogicalVolume* fWirePlane1Logical;
    G4LogicalVolume* fWirePlane2Logical;
    G4LogicalVolume* fCellLogical;
    G4LogicalVolume* fHadCalScintiLogical;
    G4LogicalVolume* fMagneticLogical;
    
    std::vector<G4VisAttributes*> fVisAttributes;
    
    G4double fArmAngle;
    G4RotationMatrix* fArmRotation;
    G4VPhysicalVolume* fSecondArmPhys;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
