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
// $Id: G4ImportanceGeometryConstructor.hh,v 1.2 2002-07-12 10:40:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ImportanceGeometryConstructor
//
// Class description:
//
// This class constructs a geometry which may be importance values
// asigned to. 
// The geometries are verry simple. The solid of the world
// geometry may be choosen by giving the name of the solid.
// Currently only "tube" is suported. A tube is constructed
// along the z-axis.
// The construction of volumes inside the world geometry
// is delegated to solid type specific construction classes.
// In the case of a tube geometry "cells" inside the world are 
// constructed with  G4ITubeFactory. Whcih is able to create
// slabs ranging between given values on the z-axis and filling
// the cross section of the world tube.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ImportanceGeometryConstructor_hh
#define G4ImportanceGeometryConstructor_hh G4ImportanceGeometryConstructor_hh

#include "globals.hh"
#include "g4std/set"

typedef G4std::set<G4String > G4ImportanceSolidTypes;

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4UImessenger;
class G4Material;
class G4VIStore;
class G4ImportanceGeometryMessenger;
class G4WorldImpMess;

class G4ImportanceGeometryConstructor {
public:
  G4ImportanceGeometryConstructor();
  ~G4ImportanceGeometryConstructor(){};

  void SetSolidType(const G4String &solidtype);
    // specify which type of solid the geometry
    // should consist of. At the moment only "tube" is allowed.
    // A G4IWorldTubeMessenger is created which allows to
    // define the world tube

  void ConstructWorldVolume(G4VSolid *);
    // give the solid to construct the physical
    // and logical volumes. used by G4IWorldTubeMessenger.

  void SetWorldImportance(G4double b, G4double e);
    // set the importance as i = pow(b,e)

  G4VIStore *GetIStore();
    // get the importance store for the volumes in that geometry

  G4VPhysicalVolume *GetWorldVolume();
  G4LogicalVolume *GetLogicalWorld();
  
private:

  void Error(const G4String &m){
    G4Exception("Error: G4ImportanceGeometryConstructor: " + m);
  }

  G4ImportanceGeometryMessenger *fGeoMessenger;

  G4String fSolidType;
  G4UImessenger *fSolidMessenger;
  G4Material *fGalactic;
  G4VSolid *fWorldSolid;
  G4VPhysicalVolume *fWorldVolume;
  G4LogicalVolume *fLogicalWorld;
  G4VIStore *fIStore;
  G4ImportanceSolidTypes fSolidTypes;  
  G4WorldImpMess *fWorldImpMess;
};

#endif
