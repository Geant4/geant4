// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02DetectorConstruction.hh,v 1.3 2000-12-04 16:24:05 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

#ifndef ExN02DetectorConstruction_h
#define ExN02DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "ExN02MagneticField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExN02DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

class ExN02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     ExN02DetectorConstruction();
    ~ExN02DetectorConstruction();

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
     
  private:
    
     G4Box*             solidWorld;    //pointer to the solid enveloppe 
     G4LogicalVolume*   logicWorld;    //pointer to the logical enveloppe
     G4VPhysicalVolume* physiWorld;    //pointer to the physical enveloppe
     
     G4Box*             solidTarget;   //pointer to the solid Target
     G4LogicalVolume*   logicTarget;   //pointer to the logical Target
     G4VPhysicalVolume* physiTarget;   //pointer to the physical Target
               
     G4Box*             solidTracker;  //pointer to the solid Tracker
     G4LogicalVolume*   logicTracker;  //pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker;  //pointer to the physical Tracker
     
     G4Box*             solidChamber;  //pointer to the solid Chamber
     G4LogicalVolume*   logicChamber;  //pointer to the logical Chamber
     G4VPhysicalVolume* physiChamber;  //pointer to the physical Chamber
     
     G4Material*         TargetMater;  //pointer to the target  material
     G4Material*         ChamberMater; //pointer to the chamber material     
     ExN02MagneticField* fpMagField;   //pointer to the magnetic field 
     
     ExN02DetectorMessenger* detectorMessenger;  //pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume
     G4double fTargetLength;           // Full length of Target
     G4double fTrackerLength;          // Full length of Tracker
     G4int    NbOfChambers;            // Nb of chambers in the tracker region
     G4double ChamberWidth;            // width of the chambers
     G4double ChamberSpacing;	       // distance between chambers
};

#endif

