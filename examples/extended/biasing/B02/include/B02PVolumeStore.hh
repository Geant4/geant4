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
// $Id: B02PVolumeStore.hh,v 1.1 2002-11-08 14:52:16 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class B02PVolumeStore
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef B02PVolumeStore_hh
#define B02PVolumeStore_hh B02PVolumeStore_hh

#include "globals.hh"
#include "g4std/set"
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

typedef G4std::set< G4GeometryCell, G4GeometryCellComp > B02SetGeometryCell;

class B02PVolumeStore {
public:
  B02PVolumeStore();
  ~B02PVolumeStore();
  
  void AddPVolume(const G4GeometryCell &cell);
  const G4VPhysicalVolume *GetPVolume(const G4String &name) const;
  G4String GetPNames() const;

private:
  B02SetGeometryCell fSetGeometryCell;
};



#endif
