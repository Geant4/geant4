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
// $Id: pymodG4geometry.cc,v 1.4 2006/06/29 15:33:00 gunter Exp $
// $Name: geant4-09-00 $
// ====================================================================
//   pymodG4geometry.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4VTouchable();
void export_G4TouchableHistory();
void export_G4VPhysicalVolume();
void export_G4PVPlacement();
void export_G4PVReplica();
void export_G4LogicalVolume();
void export_G4Region();
void export_G4VSolid();
void export_G4Box();
void export_G4Cons();
void export_G4Para();
void export_G4Torus();
void export_G4Trd();
void export_G4Orb();
void export_G4Sphere();
void export_G4Trap();
void export_G4Tubs();
void export_G4TransportationManager();
void export_G4FieldManager();
void export_G4Field();
void export_G4MagneticField();
void export_G4UniformMagField();
void export_G4ChordFinder();

BOOST_PYTHON_MODULE(G4geometry)
{
  export_G4VTouchable();
  export_G4TouchableHistory();
  export_G4VPhysicalVolume();
  export_G4PVPlacement();
  export_G4PVReplica();
  export_G4LogicalVolume();
  export_G4Region();
  export_G4VSolid();
  export_G4Box();
  export_G4Cons();
  export_G4Para();
  export_G4Torus();
  export_G4Trd();
  export_G4Orb();
  export_G4Sphere();
  export_G4Trap();
  export_G4Tubs();
  export_G4TransportationManager();
  export_G4FieldManager();
  export_G4Field();
  export_G4MagneticField();
  export_G4UniformMagField();
  export_G4ChordFinder();
}

