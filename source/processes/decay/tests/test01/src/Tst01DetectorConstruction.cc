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
// $Id: Tst01DetectorConstruction.cc,v 1.2 2001-07-11 10:02:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst01DetectorConstruction::Tst01DetectorConstruction()
{
}

Tst01DetectorConstruction::~Tst01DetectorConstruction()
{
}



G4VPhysicalVolume* Tst01DetectorConstruction::Construct()
{ 
  G4double density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* vacuum = new G4Material("Vacuum", 1., 1.01*g/mole, density,kStateGas,temperature,pressure);

  G4Box * mySimpleBox = new G4Box("SBox", 1.0*km, 1.0*km, 1.0*km);
  G4LogicalVolume * simpleBoxLog = new G4LogicalVolume( mySimpleBox,vacuum ,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),"SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

