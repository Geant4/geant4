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
// $Id: B01DetectorConstruction.cc,v 1.5 2002-09-17 13:59:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "B01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "g4std/strstream"
#include "PhysicalConstants.h"

#include "globals.hh"

B01DetectorConstruction::
B01DetectorConstruction(G4VPhysicalVolume &worldvol)
 : fWorldVolume(&worldvol)
{;}

B01DetectorConstruction::~B01DetectorConstruction()
{;}

G4VPhysicalVolume* B01DetectorConstruction::Construct()
{
  return fWorldVolume;
}
