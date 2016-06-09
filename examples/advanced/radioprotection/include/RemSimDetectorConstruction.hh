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
//    **************************************
//    *                                    *
//    *    RemSimDetectorConstruction.hh   *
//    *                                    *          
//    **************************************
//
// $Id: RemSimDetectorConstruction.hh,v 1.13 2006-06-29 16:22:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
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
  void AddShelterSPE(G4String);
  void AddHabitatRoof(G4String);
  void ChangeShieldingThickness(G4double);
  void ChangeRoofThickness(G4double);
  void ChooseConfiguration(G4String);

private:

  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* phantomPhys;
  G4VPhysicalVolume* detectorPhys; 
  RemSimVGeometryComponent* geometry;
  RemSimMaterial*  pMaterial;
  RemSimDetectorMessenger* messenger; 
  G4String decoratorValue;
  G4String astronautValue;
  RemSimDecorator* decorator; 
  RemSimDecorator* decoratorSPE;
  RemSimDecorator* decorator1;
  RemSimDecorator* decoratorRoof;
  G4bool moon;
  G4bool flag;
};
#endif

