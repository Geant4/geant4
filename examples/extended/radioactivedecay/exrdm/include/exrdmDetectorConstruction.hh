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
#ifndef exrdmDetectorConstruction_h
#define exrdmDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
//#include "exrdmMagneticField.hh"

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class exrdmDetectorMessenger;
class exrdmMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     exrdmDetectorConstruction();
    ~exrdmDetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     const 
     G4VPhysicalVolume* GetDetector() {return physiDetector;};

     G4double GetDetectoFullLength() {return fDetectorLength;};
     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     G4double GetDetectorThickness() {return fDetectorThickness;};
     G4double GetTargetRadius()  {return fTargetRadius;};
     G4double GetWorldRadius()   {return fWorldRadius;}; 
     
     void setTargetMaterial (G4String);
     void setDetectorMaterial(G4String);

     void setTargetRadius (G4double value) { fTargetRadius = value; };
     void setDetectorThickness(G4double value){ fDetectorThickness = value;};  
     void setTargetLength (G4double value) { fTargetLength = value; };
     void setDetectorLength(G4double value){ fDetectorLength = value;};  
     
  private:
     void DefineMaterials();
     
  private: 

     G4Tubs*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
     
     G4Tubs*             solidTarget;   // pointer to the solid Target
     G4LogicalVolume*   logicTarget;   // pointer to the logical Target
     G4VPhysicalVolume* physiTarget;   // pointer to the physical Target
               
     G4Tubs*             solidDetector;  // pointer to the solid Detector
     G4LogicalVolume*   logicDetector;  // pointer to the logical Detector
     G4VPhysicalVolume* physiDetector;  // pointer to the physical Detector
     
     exrdmDetectorMessenger* detectorMessenger;  // pointer to the Messenger
     exrdmMaterial* materialsManager;         // material manager
      
     G4Material* DefaultMater;          // Default material
     G4Material* TargetMater;           // Target material
     G4Material* DetectorMater;         // Detector material
 
     G4double fWorldLength;            // Full length the world volume
     G4double fTargetLength;           // Full length of the target
     G4double fDetectorLength;         // Full length of the Detector
     G4double fWorldRadius;            // Radius of  the world volume
     G4double fTargetRadius;           // Radius of the target
     G4double fDetectorThickness;      // Thickness of the Detector

     G4Region*   targetRegion;
     G4Region*   detectorRegion;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
