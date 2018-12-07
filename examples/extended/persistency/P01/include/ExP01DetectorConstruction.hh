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
/// \file persistency/P01/include/ExP01DetectorConstruction.hh
/// \brief Definition of the ExP01DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExP01DetectorConstruction_h
#define ExP01DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "ExP01MagneticField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExP01DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector Construction for the persistency example

class ExP01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     ExP01DetectorConstruction();
    ~ExP01DetectorConstruction();

  public:
  
     virtual G4VPhysicalVolume* Construct();
     
     const 
     G4VPhysicalVolume* GetTracker() {return fPhysiTracker;};
     G4double GetTrackerFullLength() {return fTrackerLength;};
     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     
     void SetTargetMaterial (G4String);
     void SetChamberMaterial(G4String);
     void SetMagField(G4double);
     
  private:

     G4Box*             fSolidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   fLogicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* fPhysiWorld;    // pointer to the physical envelope
     
     G4Box*             fSolidTarget;   // pointer to the solid Target
     G4LogicalVolume*   fLogicTarget;   // pointer to the logical Target
     G4VPhysicalVolume* fPhysiTarget;   // pointer to the physical Target
               
     G4Box*             fSolidTracker;  // pointer to the solid Tracker
     G4LogicalVolume*   fLogicTracker;  // pointer to the logical Tracker
     G4VPhysicalVolume* fPhysiTracker;  // pointer to the physical Tracker
     
     G4Box*             fSolidChamber;  // pointer to the solid Chamber
     G4LogicalVolume*   fLogicChamber;  // pointer to the logical Chamber
     G4VPhysicalVolume* fPhysiChamber;  // pointer to the physical Chamber
     
     G4Material*         fTargetMater;  // pointer to the target  material
     G4Material*         fChamberMater; // pointer to the chamber material     
     ExP01MagneticField* fPMagField;   // pointer to the magnetic field 
     
     ExP01DetectorMessenger* fDetectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume
     G4double fTargetLength;           // Full length of Target
     G4double fTrackerLength;          // Full length of Tracker
     G4int    fNbOfChambers;            // Nb of chambers in the tracker region
     G4double fChamberWidth;            // width of the chambers
     G4double fChamberSpacing;          // distance between chambers
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
