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
// $Id: Tst33IStoreBuilder.cc,v 1.15 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33IStoreBuilder.cc
//
// ----------------------------------------------------------------------

#include "Tst33IStoreBuilder.hh"
#include "Tst33CellScorerStore.hh"
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

  Tst33CellScorerStore tst33store;

  const G4VPhysicalVolume &pworld = samplegeo->GetWorldVolumeAddress();
  G4IStore *istore=0;
  istore = new G4IStore(pworld);
  if (!istore) {
    G4Exception("Tst33IStoreBuilder::CreateIStore()", "TST33-04", 
                   FatalException, " new failed to create G4IStore!");
    }
  // adding GeometryCell for world volume. ReplicaNumer = 0, since  "geomvol-V05-00-01 !
  G4GeometryCell gWorldCell(pworld, 0);
  istore->AddImportanceGeometryCell(1, gWorldCell);

  tst33store.AddG4CellScorer(gWorldCell);
  
  G4int i=1;
  for (i=1; i <= 19; ++i) {
    G4double imp = std::pow(2.0,i-1);
    G4GeometryCell gCell(samplegeo->GetGeometryCell(i, ""));
    if (i==19) {
	imp = std::pow(2.0,17);
    }
    else {
      G4GeometryCell gCellMinus(samplegeo->GetGeometryCell(i, "I1-"));
      G4GeometryCell gCellPlus(samplegeo->GetGeometryCell(i, "I1+"));
    
//       istore->AddImportanceGeometryCell(imp, gCellMinus);
//       istore->AddImportanceGeometryCell(imp, gCellPlus);
      istore->AddImportanceGeometryCell(imp, gCellMinus.GetPhysicalVolume(), i); //ASO - change overloading of AddImportanceGeometryCell?
      istore->AddImportanceGeometryCell(imp, gCellPlus.GetPhysicalVolume(), i); //ASO
      tst33store.AddG4CellScorer(gCell);
    }
    //    istore->AddImportanceGeometryCell(imp, gCell);
    istore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), i);
    G4cout << " adding importance to: " << gCell.GetPhysicalVolume().GetName() << G4endl;

  }
  return istore;
}

