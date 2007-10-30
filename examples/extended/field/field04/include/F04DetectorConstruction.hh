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
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F04DetectorConstruction_h
#define F04DetectorConstruction_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4RotationMatrix.hh"

class G4Tubs;

class G4LogicalVolume;
class G4VPhysicalVolume;

class F04Materials;
class G4Material;

class F04SimpleSolenoid;
class F04FocusSolenoid;

class F04DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"

class F04DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    F04DetectorConstruction();
    ~F04DetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    void UpdateGeometry();

    // stringToRotationMatrix() converts a string "X90,Y45" into a
    // G4RotationMatrix.
    // This is an active rotation, in that the object is first rotated
    // around the parent's X axis by 90 degrees, then the object is
    // further rotated around the parent's Y axis by 45 degrees.
    // The return value points to a G4RotationMatrix on the heap, so
    // it is persistent. Angles are in degrees, can have decimals,
    // and can be negative. Axes are X, Y, Z.

    static G4RotationMatrix stringToRotationMatrix(G4String rotation);

  public:

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);

     void SetCaptureMgntRadius(G4double);
     void SetCaptureMgntLength(G4double);
     void SetCaptureMgntB1(G4double);
     void SetCaptureMgntB2(G4double);

     void SetTransferMgntRadius(G4double);
     void SetTransferMgntLength(G4double);
     void SetTransferMgntB(G4double);
     void SetTransferMgntPos(G4double);

     void SetTargetMaterial (G4String);
     void SetTargetThickness(G4double);
     void SetTargetRadius(G4double);
     void SetTargetPos(G4double);
     void SetTargetAngle(G4int);
     
     void SetDegraderMaterial (G4String);     
     void SetDegraderThickness(G4double);     
     void SetDegraderRadius(G4double);
     void SetDegraderPos(G4double);          
      
  public:
  
     G4Material* GetWorldMaterial()    {return WorldMaterial;};
     G4double GetWorldSizeZ()          {return WorldSizeZ;}; 
     G4double GetWorldSizeR()          {return WorldSizeR;};

     G4double GetCaptureMgntRadius()   {return CaptureMgntRadius;};
     G4double GetCaptureMgntLength()   {return CaptureMgntLength;};
     G4double GetCaptureMgntB1()       {return CaptureMgntB1;};
     G4double GetCaptureMgntB2()       {return CaptureMgntB2;};

     G4double GetTransferMgntRadius()  {return TransferMgntRadius;};
     G4double GetTransferMgntLength()  {return TransferMgntLength;};
     G4double GetTransferMgntB()       {return TransferMgntB;};
     G4double GetTransferMgntPos()     {return TransferMgntPos;};

     G4Material* GetTargetMaterial()  {return TargetMaterial;};
     G4double    GetTargetRadius()    {return TargetRadius;};
     G4double    GetTargetThickness() {return TargetThickness;};
     G4double    GetTargetPos()       {return TargetPos;};
     G4int       GetTargetAngle()     {return TargetAngle;};

     G4Material* GetDegraderMaterial()  {return DegraderMaterial;};
     G4double    GetDegraderRadius()    {return DegraderRadius;};
     G4double    GetDegraderThickness() {return DegraderThickness;};
     G4double    GetDegraderPos()       {return DegraderPos;}; 

  private:

     F04Materials* materials;

     G4Material* Vacuum;
     
     G4Tubs*            solidWorld;
     G4LogicalVolume*   logicWorld;
     G4VPhysicalVolume* physiWorld;

     G4Tubs*            solidTarget;
     G4LogicalVolume*   logicTarget;
     G4VPhysicalVolume* physiTarget;

     G4Tubs*            solidDegrader;
     G4LogicalVolume*   logicDegrader;
     G4VPhysicalVolume* physiDegrader;

     G4Tubs*            solidCaptureMgnt;
     G4LogicalVolume*   logicCaptureMgnt;
     G4VPhysicalVolume* physiCaptureMgnt;

     G4Tubs*            solidTransferMgnt;
     G4LogicalVolume*   logicTransferMgnt;
     G4VPhysicalVolume* physiTransferMgnt;

     G4Material*        WorldMaterial;
     G4double           WorldSizeR;
     G4double           WorldSizeZ;

     G4double           CaptureMgntLength;
     G4double           CaptureMgntRadius;
     G4double           CaptureMgntB1;
     G4double           CaptureMgntB2;

     G4double           TransferMgntLength;
     G4double           TransferMgntRadius;
     G4double           TransferMgntB;
     G4double           TransferMgntPos;

     G4Material*        TargetMaterial;
     G4double           TargetThickness;
     G4double           TargetRadius;
     G4double           TargetPos;
     G4int              TargetAngle;

     G4Material*        DegraderMaterial;
     G4double           DegraderThickness;
     G4double           DegraderRadius;
     G4double           DegraderPos;

     F04FocusSolenoid*  focusSolenoid;
     F04SimpleSolenoid* simpleSolenoid;

  private:
    
     void DefineMaterials();

     F04DetectorMessenger* detectorMessenger;  // pointer to the Messenger

};

#endif
