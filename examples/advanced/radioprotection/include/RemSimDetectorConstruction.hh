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
// $Id: RemSimDetectorConstruction.hh,v 1.6 2004-05-21 11:14:51 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef RemSimDetectorConstruction_H
#define RemSimDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class RemSimMaterial;
class RemSimDetectorMessenger;
class G4Box;
class G4VisAttributes;
class RemSimVGeometryComponent;
class RemSimDecorator;
class RemSimDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  RemSimDetectorConstruction();
  ~RemSimDetectorConstruction();

  G4VPhysicalVolume* Construct();

  void ConstructVolume();
  void AddShielding(G4String);
  void AddAstronaut();
  void ChangeShieldingMaterial(G4String);
  void ChangeShieldingThickness(G4double);
  void ChooseConfiguration(G4String);

private:
  RemSimDecorator* retrieveDecorator();
  RemSimDecorator* retrieveDecorator1();

  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* phantomPhys;
  G4VPhysicalVolume* detectorPhys;
  RemSimVGeometryComponent* pVehicle;
  RemSimVGeometryComponent* pMoon;
  RemSimMaterial*  pMaterial;
  RemSimDetectorMessenger* messenger; 
  RemSimDecorator* decorator;
  RemSimDecorator* decorator1; 
  G4String decoratorValue;
  G4String astronautValue;
  G4bool moon;
};
#endif

