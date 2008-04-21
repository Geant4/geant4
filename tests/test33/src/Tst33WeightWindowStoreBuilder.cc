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
// $Id: Tst33WeightWindowStoreBuilder.cc,v 1.7 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33WeightWindowStoreBuilder.cc
//
// ----------------------------------------------------------------------

#include "Tst33WeightWindowStoreBuilder.hh"
#include "G4WeightWindowStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryCell.hh"
#include "globals.hh"
#include "Tst33VGeometry.hh"

Tst33WeightWindowStoreBuilder::Tst33WeightWindowStoreBuilder()
{}

Tst33WeightWindowStoreBuilder::~Tst33WeightWindowStoreBuilder()
{}

G4VWeightWindowStore *Tst33WeightWindowStoreBuilder::CreateWeightWindowStore(Tst33VGeometry *samplegeo) {
  // create an importance store and fill it with the importance
  // per cell values
  const G4VPhysicalVolume &pworld = samplegeo->GetWorldVolumeAddress();
  G4cout << " weight window store name: " << pworld.GetName() << G4endl;
  G4WeightWindowStore *wwstore=0;
  wwstore = new G4WeightWindowStore(pworld);


  // weights for the world volume, general energy bounds

  G4GeometryCell gWorldCell(pworld, 0);

  std::set<G4double, std::less<G4double> > enBounds;
  enBounds.insert(100 * MeV);
  enBounds.insert(10 * MeV);
  enBounds.insert(1 * MeV);

  wwstore->SetGeneralUpperEnergyBounds(enBounds);

  std::vector<G4double> lowerWeightsWorld;
  lowerWeightsWorld.push_back(1);
  lowerWeightsWorld.push_back(1);
  lowerWeightsWorld.push_back(1);
  
  wwstore->AddLowerWeights(gWorldCell, lowerWeightsWorld);

  
  G4int i(1);
  for (i=1; i <= 19; ++i) {
    G4double lowerWeight = 1./std::pow(2.0,i-1);
    G4GeometryCell gCell(samplegeo->GetGeometryCell(i, ""));

    std::vector<G4double> lowerWeights;

    if (i==19) {
	lowerWeight = 1./std::pow(2.0,17);
    }

    lowerWeights.push_back(lowerWeight);
    lowerWeights.push_back(lowerWeight);
    lowerWeights.push_back(lowerWeight);

    if (i!=19)  {
      G4GeometryCell gCellMinus(samplegeo->GetGeometryCell(i, "I1-"));
      G4GeometryCell gCellPlus(samplegeo->GetGeometryCell(i, "I1+"));
      wwstore->AddLowerWeights(G4GeometryCell(gCellMinus.GetPhysicalVolume(),i), lowerWeights);
      wwstore->AddLowerWeights(G4GeometryCell(gCellPlus.GetPhysicalVolume(),i), lowerWeights);
    }

    wwstore->AddLowerWeights(G4GeometryCell(gCell.GetPhysicalVolume(),i), lowerWeights);

  }
  return wwstore;
}

