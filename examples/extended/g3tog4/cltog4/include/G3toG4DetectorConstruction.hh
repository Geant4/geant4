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
// $Id: G3toG4DetectorConstruction.hh,v 1.2 2001-07-11 09:58:10 gunter Exp $
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


