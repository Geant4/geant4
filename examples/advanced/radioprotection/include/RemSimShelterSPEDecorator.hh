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
// $Id: RemSimShelterSPEDecorator.hh,v 1.2 2004/05/22 12:57:05 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#ifndef RemSimShelterSPEDecorator_h
#define RemSimShelterSPEDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
#include "globals.hh"
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;
class G4VPhysicalVolume;
class G4VisAttributes;

class RemSimShelterSPEDecorator: public RemSimDecorator
{
public:
  RemSimShelterSPEDecorator(RemSimVGeometryComponent*);
  ~RemSimShelterSPEDecorator();

  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void PrintDetectorParameters();
  G4VPhysicalVolume* GetShelter(){return 0;};
  void ChangeMother(G4VPhysicalVolume*){;};

private:
  void ConstructShelterSPE(G4VPhysicalVolume*);

  G4double shelterSPEX;
  G4double shelterSPEY;
  G4double shelterSPEZ;
  G4double translation;
  RemSimMaterial* pMaterial; 
  G4String shelterSPEMaterial;
  G4VisAttributes* shelterSPEVisAtt;
  G4Box* shelterSPE;
  G4LogicalVolume* shelterSPELog;
  G4VPhysicalVolume* shelterSPEPhys;
};
#endif
