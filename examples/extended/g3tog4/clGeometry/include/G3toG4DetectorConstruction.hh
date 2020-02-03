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
/// \file g3tog4/clGeometry/include/G3toG4DetectorConstruction.hh
/// \brief Definition of the G3toG4DetectorConstruction class
//
//
//
#ifndef G3toG4DetectorConstruction_h
#define G3toG4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G3G4Interface.hh"
#include "globals.hh"

/// Detector construction class. 
///
/// Most the work is done in G4BuildGeom(), which returns a pointer to 
/// the top-level logical volume in the detector defined by 
/// the call list file inFile.

class G3toG4DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  G3toG4DetectorConstruction(G4String inFile="svt.dat");
  virtual ~G3toG4DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
 
private:
  G4LogicalVolume* SimpleConstruct();
  G4String fInFile;
};

#endif


