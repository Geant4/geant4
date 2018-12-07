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
/// \file field/field04/include/F04DetectorConstruction.hh
/// \brief Definition of the F04DetectorConstruction class
//

#ifndef F04DetectorConstruction_h
#define F04DetectorConstruction_h 1

#include "globals.hh"
#include "G4Cache.hh"

#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4Tubs;

class G4VPhysicalVolume;

class F04Materials;
class G4Material;

class F04GlobalField;

class F04DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"

class F04DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    F04DetectorConstruction();
    virtual ~F04DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    virtual void ConstructSDandField();

    // StringToRotationMatrix() converts a string "X90,Y45" into a
    // G4RotationMatrix.
    // This is an active rotation, in that the object is first rotated
    // around the parent's X axis by 90 degrees, then the object is
    // further rotated around the parent's Y axis by 45 degrees.
    // The return value points to a G4RotationMatrix on the heap, so
    // it is persistent. Angles are in degrees, can have decimals,
    // and can be negative. Axes are X, Y, Z.

    static G4RotationMatrix StringToRotationMatrix(G4String rotation);

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
 
     G4Material* GetWorldMaterial()    {return fWorldMaterial;}
     G4double GetWorldSizeZ()          {return fWorldSizeZ;}
     G4double GetWorldSizeR()          {return fWorldSizeR;}

     G4LogicalVolume* GetCaptureMgnt()     {return fLogicCaptureMgnt;}
     G4double GetCaptureMgntRadius()       {return fCaptureMgntRadius;}
     G4double GetCaptureMgntLength()       {return fCaptureMgntLength;}
     G4double GetCaptureMgntB1()           {return fCaptureMgntB1;}
     G4double GetCaptureMgntB2()           {return fCaptureMgntB2;}
     G4ThreeVector GetCaptureMgntCenter()  {return fCaptureMgntCenter;}

     G4LogicalVolume* GetTransferMgnt()    {return fLogicTransferMgnt;}
     G4double GetTransferMgntRadius()      {return fTransferMgntRadius;}
     G4double GetTransferMgntLength()      {return fTransferMgntLength;}
     G4double GetTransferMgntB()           {return fTransferMgntB;}
     G4double GetTransferMgntPos()         {return fTransferMgntPos;}
     G4ThreeVector GetTransferMgntCenter() {return fTransferMgntCenter;}

     G4Material* GetTargetMaterial()  {return fTargetMaterial;}
     G4double    GetTargetRadius()    {return fTargetRadius;}
     G4double    GetTargetThickness() {return fTargetThickness;}
     G4double    GetTargetPos()       {return fTargetPos;}
     G4int       GetTargetAngle()     {return fTargetAngle;}

     G4Material* GetDegraderMaterial()  {return fDegraderMaterial;}
     G4double    GetDegraderRadius()    {return fDegraderRadius;}
     G4double    GetDegraderThickness() {return fDegraderThickness;}
     G4double    GetDegraderPos()       {return fDegraderPos;}

  private:

     F04DetectorMessenger* fDetectorMessenger;  // pointer to the Messenger
     G4Cache<F04GlobalField*> fFieldSetUp;

     F04Materials* fMaterials;

     G4Material* fVacuum;
 
     G4Tubs*            fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;

     G4Tubs*            fSolidTarget;
     G4LogicalVolume*   fLogicTarget;
     G4VPhysicalVolume* fPhysiTarget;

     G4Tubs*            fSolidDegrader;
     G4LogicalVolume*   fLogicDegrader;
     G4VPhysicalVolume* fPhysiDegrader;

     G4Tubs*            fSolidCaptureMgnt;
     G4LogicalVolume*   fLogicCaptureMgnt;
     G4VPhysicalVolume* fPhysiCaptureMgnt;

     G4Tubs*            fSolidTransferMgnt;
     G4LogicalVolume*   fLogicTransferMgnt;
     G4VPhysicalVolume* fPhysiTransferMgnt;

     G4Material*        fWorldMaterial;
     G4double           fWorldSizeR;
     G4double           fWorldSizeZ;

     G4double           fCaptureMgntLength;
     G4double           fCaptureMgntRadius;
     G4double           fCaptureMgntB1;
     G4double           fCaptureMgntB2;

     G4double           fTransferMgntLength;
     G4double           fTransferMgntRadius;
     G4double           fTransferMgntB;
     G4double           fTransferMgntPos;

     G4Material*        fTargetMaterial;
     G4double           fTargetThickness;
     G4double           fTargetRadius;
     G4double           fTargetPos;
     G4int              fTargetAngle;

     G4Material*        fDegraderMaterial;
     G4double           fDegraderThickness;
     G4double           fDegraderRadius;
     G4double           fDegraderPos;

     G4ThreeVector fCaptureMgntCenter, fTransferMgntCenter;

  private:

     void DefineMaterials();

};

#endif
