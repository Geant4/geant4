#ifndef BrachyDetectorConstructionLeipzig_H
#define BrachyDetectorConstructionLeipzig_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"


class G4LogicalVolume;

class G4Tubs;
class G4Box;
class G4Sphere;
class G4Tubs;
class G4Colour;
class G4VPhysicalVolume;
class G4VPhysicalVolume;
class BrachyMaterial;

class BrachyDetectorConstructionLeipzig
{
 public:
	BrachyDetectorConstructionLeipzig();
	~BrachyDetectorConstructionLeipzig();


 private:

   G4VisAttributes* simpleCapsuleVisAtt;
   G4VisAttributes* simpleCapsuleTipVisAtt;
public:
  
  void  ConstructLeipzig(G4VPhysicalVolume*);

 private:G4Tubs* Capsule ;
  G4LogicalVolume*  CapsuleLog;    //pointer to the logical World
  G4VPhysicalVolume* CapsulePhys;
 
  G4Sphere* CapsuleTip;
  G4LogicalVolume* CapsuleTipLog;
  G4VPhysicalVolume* CapsuleTipPhys;

  G4Tubs* IridiumCore;
  G4LogicalVolume* IridiumCoreLog;
  G4VPhysicalVolume* IridiumCorePhys;
  
  G4Tubs* Appl1;
  G4LogicalVolume* Appl1Log ;
  G4VPhysicalVolume* Appl1Phys;

  G4Tubs* Appl2;
  G4LogicalVolume* Appl2Log ;
  G4VPhysicalVolume* Appl2Phys;

  BrachyMaterial* pMat;   
  };

#endif
