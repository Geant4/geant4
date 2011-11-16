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
// $Id: Tst33ScorerBuilder.cc,v 1.9 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33ScorerBuilder.cc
//
// ----------------------------------------------------------------------

#include "Tst33ScorerBuilder.hh"
#include "G4CellScorerStore.hh"
#include "G4CellScorer.hh"
#include "Tst33VGeometry.hh"

//ASO
// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
//#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"


Tst33ScorerBuilder::Tst33ScorerBuilder()
{}

Tst33ScorerBuilder::~Tst33ScorerBuilder()
{}


G4CellScorerStore *Tst33ScorerBuilder::
CreateScorer(Tst33VGeometry *samplegeo, 
	     const G4CellScorer **specialCellScorer){


  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //
  // Sensitive Detector Name
  G4String phantomSDname = "PhantomSDColl";
  G4String phantomSDnamePlus = "PhantomSDPlus";
  G4String phantomSDnameMinus = "PhantomSDMinus";

  //------------------------
  // MultiFunctionalDetector
  //------------------------
  //
  // Define MultiFunctionalDetector with name.
  G4MultiFunctionalDetector* MFDet2 = new G4MultiFunctionalDetector(phantomSDname);
  G4MultiFunctionalDetector* MFDetPlus = new G4MultiFunctionalDetector(phantomSDnamePlus);
  G4MultiFunctionalDetector* MFDetMinus = new G4MultiFunctionalDetector(phantomSDnameMinus);
  SDman->AddNewDetector( MFDet2 );                 // Register SD to SDManager
  SDman->AddNewDetector( MFDetPlus );                 // Register SD to SDManager
  SDman->AddNewDetector( MFDetMinus );                 // Register SD to SDManager

  G4GeometryCell gWorldCell(samplegeo->GetWorldVolumeAddress(), -1);
  
  G4CellScorerStore *cs_store = new G4CellScorerStore();
  if (!cs_store) {
    G4Exception("Tst33ScorerBuilder::CreateScorer()",
        "TST33-08", FatalException, " new failed to create G4CellScorerStore!");
  }
  cs_store->AddCellScorer(gWorldCell);
  
  
  G4int i = 1;
  for (i=1; i <= 19; i++) {
    G4GeometryCell gCell(samplegeo->GetGeometryCell(i, ""));
    //xtest    gCell.GetPhysicalVolume().GetLogicalVolume()->SetSensitiveDetector(MFDet2);

    const G4CellScorer *s = cs_store->AddCellScorer(gCell);

    if (i==18) {
      *specialCellScorer = s; // ????
    }
    if (i!=19) {
      G4GeometryCell gCellMinus(samplegeo->GetGeometryCell(i, "I1-"));
      gCellMinus.GetPhysicalVolume().GetLogicalVolume()->SetSensitiveDetector(MFDetMinus);
      //      cs_store->AddCellScorer(gCellMinus);
      G4GeometryCell gCellPlus(samplegeo->GetGeometryCell(i, "I1+"));
      gCellPlus.GetPhysicalVolume().GetLogicalVolume()->SetSensitiveDetector(MFDetPlus);
      //      cs_store->AddCellScorer(gCellPlus);
    }
  }


  G4String psName;
  G4PSNofCollision*   scorer0coll = new G4PSNofCollision(psName="Collisions");  
  MFDet2->RegisterPrimitive(scorer0coll);

  G4PSNofCollision*   scorer0minus = new G4PSNofCollision(psName="Collisions");  
  MFDetMinus->RegisterPrimitive(scorer0minus);

  G4PSNofCollision*   scorer0plus = new G4PSNofCollision(psName="Collisions");  
  MFDetPlus->RegisterPrimitive(scorer0plus);



  G4PSNofCollision*   scorer1coll = new G4PSNofCollision(psName="CollWeight");  
  scorer1coll->Weighted(true);
  MFDet2->RegisterPrimitive(scorer1coll);

  G4PSNofCollision*   scorer1minus = new G4PSNofCollision(psName="CollWeight");  
  scorer1minus->Weighted(true);
  MFDetMinus->RegisterPrimitive(scorer1minus);

  G4PSNofCollision*   scorer1plus = new G4PSNofCollision(psName="CollWeight");  
  scorer1plus->Weighted(true);
  MFDetPlus->RegisterPrimitive(scorer1plus);


//tooAdd   G4PSPopulation*   scorer2coll = new G4PSPopulation(psName="Population");  
//tooAdd   MFDet->RegisterPrimitive(scorer2coll);

//tooAdd   G4PSPopulation*   scorer2minus = new G4PSPopulation(psName="Population");  
//tooAdd   MFDetMinus->RegisterPrimitive(scorer2minus);

//tooAdd   G4PSPopulation*   scorer2plus = new G4PSPopulation(psName="Population");  
//tooAdd   MFDetPlus->RegisterPrimitive(scorer2plus);

  G4PSTrackCounter* scorer3coll = new G4PSTrackCounter(psName="TrackEnter",fCurrent_In);  
  MFDet2->RegisterPrimitive(scorer3coll);

  G4PSTrackCounter* scorer3minus = new G4PSTrackCounter(psName="TrackEnter",fCurrent_In);  
  MFDetMinus->RegisterPrimitive(scorer3minus);

  G4PSTrackCounter* scorer3plus = new G4PSTrackCounter(psName="TrackEnter",fCurrent_In);  
  MFDetPlus->RegisterPrimitive(scorer3plus);

  G4PSTrackLength* scorer4coll = new G4PSTrackLength(psName="SL");  
  MFDet2->RegisterPrimitive(scorer4coll);
  G4PSTrackLength* scorer4minus = new G4PSTrackLength(psName="SL");  
  MFDetMinus->RegisterPrimitive(scorer4minus);
  G4PSTrackLength* scorer4plus = new G4PSTrackLength(psName="SL");  
  MFDetPlus->RegisterPrimitive(scorer4plus);

  G4PSTrackLength* scorer5coll = new G4PSTrackLength(psName="SLW");  
  scorer5coll->Weighted(true);
  MFDet2->RegisterPrimitive(scorer5coll);
  G4PSTrackLength* scorer5minus = new G4PSTrackLength(psName="SLW");  
  scorer5minus->Weighted(true);
  MFDetMinus->RegisterPrimitive(scorer5minus);
  G4PSTrackLength* scorer5plus = new G4PSTrackLength(psName="SLW");  
  scorer5plus->Weighted(true);
  MFDetPlus->RegisterPrimitive(scorer5plus);

  G4PSTrackLength* scorer6coll = new G4PSTrackLength(psName="SLWE");  
  scorer6coll->Weighted(true);
  scorer6coll->MultiplyKineticEnergy(true);
  MFDet2->RegisterPrimitive(scorer6coll);

  G4PSTrackLength* scorer6minus = new G4PSTrackLength(psName="SLWE");  
  scorer6minus->Weighted(true);
  scorer6minus->MultiplyKineticEnergy(true);
  MFDetMinus->RegisterPrimitive(scorer6minus);

  G4PSTrackLength* scorer6plus = new G4PSTrackLength(psName="SLWE");  
  scorer6plus->Weighted(true);
  scorer6plus->MultiplyKineticEnergy(true);
  MFDetPlus->RegisterPrimitive(scorer6plus);

  G4PSTrackLength* scorer7coll = new G4PSTrackLength(psName="SLW_V");  
  scorer7coll->Weighted(true);
  scorer7coll->DivideByVelocity(true);
  MFDet2->RegisterPrimitive(scorer7coll);

  G4PSTrackLength* scorer7minus = new G4PSTrackLength(psName="SLW_V");  
  scorer7minus->Weighted(true);
  scorer7minus->DivideByVelocity(true);
  MFDetMinus->RegisterPrimitive(scorer7minus);

  G4PSTrackLength* scorer7plus = new G4PSTrackLength(psName="SLW_V");  
  scorer7plus->Weighted(true);
  scorer7plus->DivideByVelocity(true);
  MFDetPlus->RegisterPrimitive(scorer7plus);

  G4PSTrackLength* scorer8coll = new G4PSTrackLength(psName="SLWE_V");  
  scorer8coll->Weighted(true);
  scorer8coll->MultiplyKineticEnergy(true);
  scorer8coll->DivideByVelocity(true);
  MFDet2->RegisterPrimitive(scorer8coll);

  G4PSTrackLength* scorer8minus = new G4PSTrackLength(psName="SLWE_V");  
  scorer8minus->Weighted(true);
  scorer8minus->MultiplyKineticEnergy(true);
  scorer8minus->DivideByVelocity(true);
  MFDetMinus->RegisterPrimitive(scorer8minus);

  G4PSTrackLength* scorer8plus = new G4PSTrackLength(psName="SLWE_V");  
  scorer8plus->Weighted(true);
  scorer8plus->MultiplyKineticEnergy(true);
  scorer8plus->DivideByVelocity(true);
  MFDetPlus->RegisterPrimitive(scorer8plus);



  return cs_store;
}
