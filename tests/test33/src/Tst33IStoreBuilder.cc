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
// $Id: Tst33IStoreBuilder.cc,v 1.3 2002-10-31 08:32:44 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33IStoreBuilder.cc
//
// ----------------------------------------------------------------------

#include "Tst33IStoreBuilder.hh"
#include "G4IStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryCell.hh"
#include "globals.hh"
#include "Tst33VGeometry.hh"

Tst33IStoreBuilder::Tst33IStoreBuilder()
{}

Tst33IStoreBuilder::~Tst33IStoreBuilder()
{}

G4VIStore *Tst33IStoreBuilder::CreateIStore(Tst33VGeometry *samplegeo) {
  // create an importance store and fill it with the importance
  // per cell values
  const G4VPhysicalVolume &pworld = samplegeo->GetWorldVolume();
  G4IStore *istore(0);
  istore = new G4IStore(pworld);
  if (!istore) {
    G4std::G4Exception("Tst33IStoreBuilder::CreateIStore new failed to create G4IStore!");
    }
  // adding GeometryCell for world volume. ReplicaNumer = -1 !
  G4GeometryCell gWorldCell(pworld, -1);
  istore->AddImportanceGeometryCell(1, gWorldCell);
  
  G4int i(1);
  for (i=1; i <= 19; ++i) {
    G4String volname = samplegeo->GetCellName(i);
    G4double imp = G4std::pow(2,i-1);
    if (i==19) {
	imp = G4std::pow(2,17);
    }
    const G4VPhysicalVolume *pvol(0);
    pvol = samplegeo->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      istore->AddImportanceGeometryCell(imp, gCell);
    }
  }
  return istore;
}

