// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02DetectorConstruction.hh,v 1.2 1999-12-15 14:49:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef ExN02DetectorConstruction_h
#define ExN02DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
// class G4UniformMagField;
#include "ExN02MagneticField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExN02DetectorMessenger;

class ExN02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN02DetectorConstruction();
    ~ExN02DetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     const G4VPhysicalVolume* GetTracker() {return physiTracker;};
     G4double GetTargetFullLength() {return fTargetLength;};
     G4double GetWorldFullLength()  {return fWorldLength;}; 
     
     void SetDetectorLength(G4double length);
     // void setMaterial(G4String val);
     void SetMagField(G4double val);
  private:
    
     G4Box*             solidWorld;    //pointer to the solid enveloppe 
     G4LogicalVolume*   logicWorld;    //pointer to the logical enveloppe
     G4VPhysicalVolume* physiWorld;    //pointer to the physical enveloppe
     G4double fullSizeWorld;           //full size of the enveloppe
          
     G4Box*             solidTracker;  //pointer to the solid Tracker
     G4LogicalVolume*   logicTracker;  //pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker;  //pointer to the physical Tracker
     
     G4Box*             solidTarget;   //pointer to the solid Target
     G4LogicalVolume*   logicTarget;   //pointer to the logical Target
     G4VPhysicalVolume* physiTarget;   //pointer to the physical Target
     
     G4Box*             solidChamber;  //pointer to the solid Chamber
     G4LogicalVolume*   logicChamber;  //pointer to the logical Chamber
     G4VPhysicalVolume* physiChamber;  //pointer to the physical Chamber
     
     G4Material*         myMaterial;   //pointer to the material
     ExN02MagneticField* fpMagField;   //pointer to the mag field 
     
     ExN02DetectorMessenger* detectorMessenger;  //pointer to the Messenger object  
     G4double fWorldLength;               // Full length of the world volume
     //  Sizes of pieces of the detector
     G4double fDetectorLength;           // Full length of the detector

     G4double fTrackerLength;            // Full length of Tracker
     G4double fTargetLength;             // Full length of Target	  

};

#endif

