//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************


#ifndef BrachyDetectorConstructionIr_H
#define BrachyDetectorConstructionIr_H 1

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

class BrachyDetectorConstructionIr 
{
 public:
	BrachyDetectorConstructionIr();
	~BrachyDetectorConstructionIr();


 private:
     
 G4VisAttributes* simpleCapsuleVisAtt;
  G4VisAttributes* simpleCapsuleTipVisAtt;
  G4VisAttributes* simpleIridiumVisAtt;
  
 

public:
  
  void  CleanIridium();
  void  ConstructIridium(G4VPhysicalVolume*);
 private:

  G4Box*              WaterBox;
  G4LogicalVolume*    WaterBoxLog;
  G4VPhysicalVolume*   WaterBoxPhys;

  G4Tubs* Capsule ;
  G4LogicalVolume*  CapsuleLog;    //pointer to the logical World
  G4VPhysicalVolume* CapsulePhys;
 
  G4Sphere* CapsuleTip;
  G4LogicalVolume* CapsuleTipLog;
  G4VPhysicalVolume* CapsuleTipPhys;

  G4Tubs* IridiumCore;
  G4LogicalVolume* IridiumCoreLog;
  G4VPhysicalVolume* IridiumCorePhys;
  
  
  BrachyMaterial* pMat;   
  };



#endif








