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
// $Id: F07DetectorConstruction.hh,v 1.10 2008-09-22 16:41:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F07DetectorConstruction_h
#define F07DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "F07MagneticField.hh"

class G4VSolid;
class G4Box;
class G4Tubs;
// class G4Cons; 
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
// class G4VPVParameterisation;
class G4UserLimits;
class F07DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F07DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     F07DetectorConstruction();
    ~F07DetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     const 
     G4VPhysicalVolume* GetTracker() {return physiTracker;};
     G4double GetTrackerFullLength() {return fTrackerLength;};
     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     
     void setTargetMaterial (G4String);
     void setChamberMaterial(G4String);
     void SetMagField(G4double);
     void SetMaxStep (G4double);     
     
  private:

     G4Box*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
     
     G4VSolid*          solidTarget;   // pointer to the solid Target
     G4LogicalVolume*   logicTarget;   // pointer to the logical Target
     G4VPhysicalVolume* physiTarget;   // pointer to the physical Target
               
     G4Tubs*            solidTracker;  // pointer to the solid Tracker
     G4LogicalVolume*   logicTracker;  // pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker;  // pointer to the physical Tracker

#if 0
     std::vector<G4Box*>             solidChambers;    // Store for the solids Chamber
     std::vector<G4VPhysicalVolume*> physicalChambers; // Store for the physical Chamber
     std::vector<G4LogicalVolume*>   logicalChambers;  // Store to the logical Chamber
#endif
  
     // Pointers to the key materials in the setup - can be changed
     G4Material*        fTargetMaterial;  
     G4Material*        fChamberAbsorberMat;
     G4Material*        fChamberGasMat;

//   G4VPVParameterisation* chamberParam; 
     G4UserLimits* stepLimit;             // pointer to user step limits

     F07MagneticField* fpMagField;   // pointer to the magnetic field 
     
     F07DetectorMessenger* detectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume
     G4double fWorldHalfWidth;         // Half width  of the world volume

     G4double fTrackerLength;          // Full Length of Tracker
     G4double fTrackerRadius;          // Full Radius of Tracker

     unsigned int fNbOfChambers;     // Nb of chambers in the tracker region
     G4double fDeltaRadius;          // Delta of radius of chambers
     G4double fForwardSphereRadius;    // Outer radius of forward chambers
     // G4double fChamberSpacing;	   // ??? distance between chambers

     G4double fTargetLength;           // Full length           of Target
     G4double fTargetMaxRadius;        // Maximum Radial extent of Target
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
