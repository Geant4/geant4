//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef test31DetectorConstruction_h
#define test31DetectorConstruction_h 1

// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31DetectorConstruction -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"
#include "globals.hh"
#include "G4ios.hh"

class test31EventAction;
class test31DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class test31DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    test31DetectorConstruction();
   ~test31DetectorConstruction();

  public:
     
    void SetAbsorberMaterial(const G4String& v) {nameMatAbsorber = v;};
    G4Material* GetAbsorberMaterial()  const {return AbsorberMaterial;};
    void SetNumberOfAbsorbers(G4int);     
    G4int GetNumberOfAbsorbers() const {return NumberOfAbsorbers;};
    void SetAbsorberThickness(G4double);     
    G4double  GetAbsorberThickness() const {return AbsorberThickness;};
    void SetAbsorberSizeXY   (G4double);            
    G4double    GetAbsorberSizeXY() const {return SizeXY;};
    void SetWorldMaterial(const G4String& v) {nameMatWorld = v;};
    G4Material* GetWorldMaterial() const {return WorldMaterial;};
    void SetWorldSizeZ   (G4double);
    G4double    GetWorldSizeZ() const {return WorldSizeZ;}; 
    void SetGap(G4double val);
    G4double    GetGap() const {return gap;}; 
    void SetNumberOfEvents(G4int val) {nEvents = val;};
    G4int GetNumberOfEvents() const {return nEvents;};
    void SetEventAction(test31EventAction* p) {theEvent = p;};
    test31EventAction* GetEventAction() const {return theEvent;};
    void SetVerbose(G4int val) {myVerbose = val;};
    G4int GetVerbose() const {return myVerbose;};
    void SetNumAbsorbersSaved(G4int val) {nAbsSaved = val;};
    G4int GetNumAbsorbersSaved() const {return nAbsSaved;};
    void SetFirstEventToDebug(G4int val) {nFirstEvtToDebug = val;};
    G4int GetFirstEventToDebug() const {return nFirstEvtToDebug;};
    void SetLastEventToDebug(G4int val) {nLastEvtToDebug = val;};
    G4int GetLastEventToDebug() const {return nLastEvtToDebug;};
    G4double  GetMaxDeltaEnergy() const {return maxDelta;};
    void SetMaxDeltaEnergy(G4double val) {maxDelta = val;};            

    const G4VPhysicalVolume* GetPhysWorld() const {return physWorld;};
    const G4LogicalVolume*   GetAbsorber()  const {return logicAbs;};

    void SetMagField(G4double,G4int);
     
    G4VPhysicalVolume* Construct();
    void UpdateGeometry();
    void PrintGeomParameters(); 
                                     
  private:

    // Methods

     void DefineMaterials();
     G4Material* GetMaterial(const G4String&);
     void ComputeGeomParameters();
     G4VPhysicalVolume* ConstructGeometry();     
     void GeometryIsChanged();     
     void MaterialIsChanged();     

    // Members
     
     G4String           nameMatAbsorber;
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;  // 
     G4double           SizeXY;
     G4double           gap;

     G4int              NumberOfAbsorbers;

     G4double           ZposAbs;
     
     G4String           nameMatWorld;
     G4Material*        WorldMaterial;
     G4double           WorldSizeZ;     
                
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physWorld;     //pointer to the physical World

     G4Box*             solidAbs;      //pointer to the solid Absorber
     G4LogicalVolume*   logicAbs;      //pointer to the logical Absorber
     G4VPhysicalVolume* physAbs;       //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     test31DetectorMessenger* detectorMessenger;   
     test31EventAction* theEvent;

     G4int myVerbose;      
     G4int nEvents;
     G4bool detIsConstructed;
     G4int nAbsSaved;
     G4int nFirstEvtToDebug;
     G4int nLastEvtToDebug;
     G4double maxDelta;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

