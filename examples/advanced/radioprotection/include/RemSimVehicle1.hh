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
// $Id: RemSimVehicle1.hh,v 1.7 2004/05/22 12:57:05 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//

#ifndef RemSimVehicle1_h
#define RemSimVehicle1_h 1

class RemSimVGeometryComponent;
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class G4VisAttributes;
class RemSimVehicle1: public RemSimVGeometryComponent
{
public:
  RemSimVehicle1();
  ~RemSimVehicle1();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  G4VPhysicalVolume* GetShelter(){return 0;};

private:
  RemSimMaterial* pMaterial;
  G4VPhysicalVolume* layervacuumPhys;  
  G4VPhysicalVolume* layer1Phys;
  G4VPhysicalVolume* layer2Phys; 
  G4VPhysicalVolume* layer3Phys;
  G4VPhysicalVolume* layer4Phys;
  G4VPhysicalVolume* layer5Phys;
  G4VPhysicalVolume* layer6Phys;
  G4VPhysicalVolume* layer7Phys;
  G4VPhysicalVolume* layer8Phys;
  G4VPhysicalVolume* layer9Phys;
  G4VPhysicalVolume* layer10Phys;
  G4VPhysicalVolume* layer11Phys; 
  G4VPhysicalVolume* layer12Phys;
  G4VPhysicalVolume* layer13Phys; 
  G4VPhysicalVolume* layer14Phys;
  G4VPhysicalVolume* layer15Phys;
  G4VPhysicalVolume* layerPhys;
  G4VPhysicalVolume* layer16Phys; 
  G4VPhysicalVolume* layer17Phys;
  G4VPhysicalVolume* layer18Phys;
  G4VPhysicalVolume* layer19Phys; 
  G4VPhysicalVolume* layer20Phys;
  G4VPhysicalVolume* layer21Phys;
  G4VPhysicalVolume* layer22Phys;
  G4VPhysicalVolume* layer23Phys;
  G4VPhysicalVolume* layer24Phys;
  G4VPhysicalVolume* layer25Phys;
  G4VPhysicalVolume* layer26Phys;
  G4VPhysicalVolume* layer27Phys;
};
#endif
