//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionI.hh     *
//    *                                      *
//    ****************************************
#ifndef BrachyDetectorConstructionI_H
#define BrachyDetectorConstructionI_H 1

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

class BrachyDetectorConstructionI
{
public:
    BrachyDetectorConstructionI();
    ~BrachyDetectorConstructionI();


private:
G4VisAttributes*   simpleMarkerVisAtt;
G4VisAttributes* simpleIodiumVisAtt;

G4VisAttributes* simpleCapsuleVisAtt;
G4VisAttributes* simpleCapsuleTipVisAtt;

public:


void  ConstructIodium(G4VPhysicalVolume*);
private:
G4Tubs* DefaultTub;
G4LogicalVolume* DefaultTub_log;
G4VPhysicalVolume*DefaultTub_Phys;

G4Tubs* Capsule ;
G4LogicalVolume*  CapsuleLog;    //pointer to the logical World
G4VPhysicalVolume* CapsulePhys;

G4Sphere* CapsuleTip;
G4LogicalVolume* CapsuleTipLog;
G4VPhysicalVolume* CapsuleTipPhys1;
G4VPhysicalVolume* CapsuleTipPhys2;

G4Tubs* IodiumCore;
G4LogicalVolume* ICoreLog;
G4VPhysicalVolume* ICorePhys;

G4Tubs* Marker;
G4LogicalVolume* MarkerLog;
G4VPhysicalVolume* MarkerPhys;

BrachyMaterial* pMat;   
};



#endif








