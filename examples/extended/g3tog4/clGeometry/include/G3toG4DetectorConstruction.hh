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
//
// $Id: G3toG4DetectorConstruction.hh,v 1.4 2006-06-29 17:20:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G3toG4DetectorConstruction_h
#define G3toG4DetectorConstruction_h 1

//--------------------------------------------------------------------------
// G3toG4DetectorConstruction. Most the work is Done in
// G4BuildGeom, which returns a G4LogicalVolume*, a pointer to the
// top-level logiical volume in the detector defined by the call List file
// inFile
//--------------------------------------------------------------------------

#include "G4VUserDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G3G4Interface.hh"
#include "globals.hh"

class G3toG4DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  G3toG4DetectorConstruction(G4String inFile="svt.dat");

  ~G3toG4DetectorConstruction();

  G4VPhysicalVolume* Construct();
  G4LogicalVolume* SimpleConstruct();

private:
  G4String _inFile;
  G4VPhysicalVolume* _pv;
  G4LogicalVolume* _lv;
};

#endif


