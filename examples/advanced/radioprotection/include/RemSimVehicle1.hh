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
// $Id$
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
