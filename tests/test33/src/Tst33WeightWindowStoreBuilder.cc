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
// $Id: Tst33WeightWindowStoreBuilder.cc,v 1.3 2003-11-25 10:20:25 gcosmo Exp $
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
  const G4VPhysicalVolume &pworld = samplegeo->GetWorldVolume();
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
    G4double lowerWeight = 1./pow(2.0,i-1);
    G4GeometryCell gCell(samplegeo->GetGeometryCell(i, ""));

    std::vector<G4double> lowerWeights;

    if (i==19) {
	lowerWeight = 1./pow(2.0,17);
    }

    lowerWeights.push_back(lowerWeight);
    lowerWeights.push_back(lowerWeight);
    lowerWeights.push_back(lowerWeight);

    if (i!=19)  {
      G4GeometryCell gCellMinus(samplegeo->GetGeometryCell(i, "I1-"));
      G4GeometryCell gCellPlus(samplegeo->GetGeometryCell(i, "I1+"));
    
      wwstore->AddLowerWeights(gCellMinus, lowerWeights);
      wwstore->AddLowerWeights(gCellPlus, lowerWeights);
    }

    wwstore->AddLowerWeights(gCell, lowerWeights);

  }
  return wwstore;
}

