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
// $Id: RemSimDetectorConstruction.hh,v 1.1 2004-01-30 12:18:24 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef RemSimDetectorConstruction_H
#define RemSimDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class RemSimVGeometryComponent;
class RemSimMaterial;
class RemSimDecorator;
class RemSimDetectorMessenger;
class RemSimDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  RemSimDetectorConstruction();
  ~RemSimDetectorConstruction();

  G4VPhysicalVolume* Construct();
  void SwitchVehicle(G4String);
  void ConstructVolume();
private:
    
  // Logical volumes
  //
  G4LogicalVolume* experimentalHall_log;
   
  // Physical volumes
  //
  G4VPhysicalVolume* experimentalHall_phys;
  G4int detectorChoice ;

  RemSimMaterial*  pMaterial;
  RemSimVGeometryComponent* pVehicle;
  RemSimDecorator* decorator;
  RemSimDetectorMessenger* messenger; 
};
#endif

