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
// $Id: tiara.cc,v 1.4 2006/06/29 15:43:16 gunter Exp $
// GEANT4 tag 
//
// 
// --------------------------------------------------------------
//      GEANT 4 - tiara
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "TiaraSim.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCutsTable.hh"
#include "TiaraTally.hh"
#include "Randomize.hh"

int main(int, char **)
{  

  std::vector<G4double> binEdges;
  binEdges.push_back(1.);
  binEdges.push_back(3.);
  binEdges.push_back(7.);
  binEdges.push_back(13.);
  binEdges.push_back(21.);
  
  TiaraTally t;
  t.setBinEdges(binEdges);
  for (G4int e=0; e<1000; e++) {    
    for (G4int i=0; i<1000; i++) {
      G4double r(1 + G4UniformRand()*20);
      t.fill(r, G4UniformRand());
    }
    t.EndOfEventAction();
  }
  for (G4int i=0; i<t.size();i++) {
    TiaraMeasure m(t.measure(i));
    G4cout << "entries: " << m.GetEntries() << G4endl;
    G4cout << m.GetMean() << G4endl;
    G4cout << m.GetVariance() << G4endl;
    G4cout << m.GetSum() << G4endl;
    G4cout << m.GetSumSquared() << G4endl;
  }

  /*
  Physics_MY_CASCADE_HP *phys = new Physics_MY_CASCADE_HP;
  phys->Construct();
  phys->SetCuts();
  G4RegionStore::GetInstance()->UpdateMaterialList();
  // Let G4ProductionCutsTable update couples
  G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable();
  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();

  if(G4ProductionCutsTable::GetProductionCutsTable()->IsModified())
  {
    phys->BuildPhysicsTable();
    G4ProductionCutsTable::GetProductionCutsTable()->PhysicsTableUpdated();
  }
  phys->DumpCutValuesTableIfRequested();
  TiaraSim &tSim = TiaraSim::GetTiaraSim();
  */
}

