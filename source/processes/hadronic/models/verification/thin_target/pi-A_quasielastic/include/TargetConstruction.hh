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
//
// $Id: TargetConstruction.hh,v 1.1 2003-07-31 01:16:41 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef TargetConstruction_h
#define TargetConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class TargetMessenger;


class TargetConstruction : public G4VUserDetectorConstruction
{
  public:
  
    TargetConstruction();
    ~TargetConstruction();

    G4VPhysicalVolume* Construct();

    void SetTargetMaterial(G4String);     
    void SetTargetThickness(G4double);     
    void SetTargetRadius(G4double);     
    void SetNumRequestedEvents(G4int);     

    void SetMagField(G4double);
     
    void UpdateGeometry();
  
    void PrintTargParameters(); 
                    
    G4Material* GetTargetMaterial()  {return TargetMaterial;};
    G4double    GetTargetThickness() {return TargetThickness;};      
    G4double    GetTargetRadius() {return TargetRadius;};      
    G4int       GetNumRequestedEvents() {return NumRequestedEvents;};      
    const G4VPhysicalVolume* GetTarget()   {return physTarg;};
     
    const G4VPhysicalVolume* GetphysWorld() {return physWorld;};           
    const G4double GetWorldLength() {return WorldLength;};           
                 
  private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructTarget();
     
  private:
     
     G4Material* WorldMaterial;
     G4Tubs* solidWorld;             // pointer to the solid world 
     G4LogicalVolume* logicWorld;    // pointer to the logical world
     G4VPhysicalVolume* physWorld;   // pointer to the physical world

     G4Material* TargetMaterial;
     G4Tubs* solidTarg;              // pointer to the solid target
     G4LogicalVolume* logicTarg;     // pointer to the logical target
     G4VPhysicalVolume* physTarg;    // pointer to the physical target
     
     G4double TargetThickness;
     G4double TargetRadius;
     G4int NumRequestedEvents;

     G4double WorldLength;

     TargetMessenger* tgtMessenger;  // pointer to the Messenger
      
};

#endif






