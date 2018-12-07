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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef eRositaDetectorConstruction_h
#define eRositaDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"


class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class eRositaDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  eRositaDetectorConstruction();
  ~eRositaDetectorConstruction();

public:
  
  G4VPhysicalVolume* Construct();
     
  const 
  G4VPhysicalVolume* GetTracker() {return physiTracker;};
  //  G4double GetTrackerFullLength() {return fTrackerLength;};
  //  G4double GetTargetFullLength()  {return fTargetLength;};
  //  G4double GetWorldFullLength()   {return fWorldLength;}; 
  
  void setTargetMaterial (G4String);
  void setTrackerMaterial(G4String);
  void setWorldMaterial(G4String);
     
private:

  G4Box*             solidWorld;    // pointer to the solid envelope 
  G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
  G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
  G4VisAttributes*   visWorld;      // pointer to visualization attributes
  
  G4Box*             solidTarget;   // pointer to the solid Target
  G4LogicalVolume*   logicTarget;   // pointer to the logical Target
  G4VPhysicalVolume* physiTarget;   // pointer to the physical Target
  G4VisAttributes*   visTarget;     // pointer to visualization attributes
  
  G4Box*             solidTracker;  // pointer to the solid Tracker
  G4LogicalVolume*   logicTracker;  // pointer to the logical Tracker
  G4VPhysicalVolume* physiTracker;  // pointer to the physical Tracker
  G4VisAttributes*   visTracker;    // pointer to visualization attributes
  
  G4Material*         TargetMater;  // pointer to the target  material
  G4Material*         TrackerMater; // pointer to the tracker material
  G4Material*         WorldMater;   // pointer to the tracker material

  G4Material*         vacuum;
 
       
  G4double hWorldLength;            // half length of the world volume
  G4double hTargetLength;           // half length of target
  G4double hTargetDepth;            // half depth of target
  G4double hTrackerLength;          // half length of tracker
  G4double hTrackerDepth;           // half depth of tracker

  G4double xPosTarget;       // x coordinate of target position
  G4double yPosTarget;       // y          -"-
  G4double zPosTarget;       // z          -"-
  G4double xPosTracker;      // x coordinate of tracker position
  G4double yPosTracker;      // y          -"-
  G4double zPosTracker;      // z          -"-
  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
