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
/// \file biasing/B03/include/B03PVolumeStore.hh
/// \brief Definition of the B03PVolumeStore class
//
//
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class B03PVolumeStore
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef B03PVolumeStore_hh
#define B03PVolumeStore_hh B03PVolumeStore_hh

#include "globals.hh"
#include <set>
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

typedef std::set< G4GeometryCell, G4GeometryCellComp > B03SetGeometryCell;

class B03PVolumeStore {
public:
  B03PVolumeStore();
  ~B03PVolumeStore();
  
  void AddPVolume(const G4GeometryCell &cell);
  const G4VPhysicalVolume *GetPVolume(const G4String &name) const;
  G4String GetPNames() const;

private:
  B03SetGeometryCell fSetGeometryCell;
};

#endif
