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
// $Id: RunAction.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation, based on 
//                         BeamTestRunAction by T. Aso
//
#include "RunAction.hh"

#include "ConfigData.hh"
#include "G4ProductionCutsTable.hh"
#include "G4RunManager.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"
#include "Run.hh"
#include <fstream>

RunAction::RunAction() 
{
  timer = new G4Timer();
}

RunAction::~RunAction() 
{
  delete timer;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  timer->Start();
}

G4Run*  RunAction::GenerateRun(){
  return new Run();
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();

  G4int numberOfReplicas = (ConfigData::GetScorerMaxTheta() - ConfigData::GetScorerMinTheta())/ConfigData::GetScorerDeltaTheta();
  
  assert (numberOfReplicas > 0);
  
  const Run* theRun = dynamic_cast<const Run*>(aRun);
  assert (0 != theRun);
  
  G4String outFile(ConfigData::GetOutputDirectory()+"/"+ConfigData::GetOutputId()+"_Result.gntplot");
  G4cout<<"Writing results to "<<outFile<<G4endl;

  std::ofstream out(outFile.data());

  std::ofstream& log = ConfigData::GetConfigFile();

  for (G4int i=0; i<theRun->GetNumberOfHitsMap(); ++i) {
    G4THitsMap<G4double>* runMap =theRun->GetHitsMap(i);
    assert (0  != runMap);
    
    if (ConfigData::GetVerbose()) {
      G4cout << " PrimitiveScorer RUN " << runMap->GetSDname() <<","<< runMap->GetName() << G4endl;
      G4cout << " Number of entries " << runMap->entries() << G4endl;
      
      std::map<G4int,G4double*>::iterator iter;
      
      for(iter = runMap->GetMap()->begin(); iter != runMap->GetMap()->end(); ++iter) {
	G4cout << "  copy no.: " << iter->first
	       << "  Run Value : " << *(iter->second)
	       << G4endl;
      }
    }
    for (G4int j=0; j<numberOfReplicas; ++j){
      G4double* flux = (*runMap)[j];
      if( !flux ) flux = new G4double(0.0);
      out << std::setw(5) << j << "\t" << std::setw(15)<< *flux*mm*mm << " mm-2 " <<G4endl;
    }
  }
  out.close();
  
  log << "===================================================================" << G4endl;
  G4ProductionCutsTable* table = G4ProductionCutsTable::GetProductionCutsTable();
  
  for (unsigned i=0; i<table->GetTableSize(); i++) {
    const G4MaterialCutsCouple* aCouple = table->GetMaterialCutsCouple(i);
    const G4ProductionCuts* aCut = aCouple->GetProductionCuts();
    
    log << " Material : " << aCouple->GetMaterial()->GetName() << G4endl;
    log << " Range cuts        : "
	<< " gamma " << G4BestUnit(aCut->GetProductionCut("gamma"),"Length")
	<< "    e- " << G4BestUnit(aCut->GetProductionCut("e-"),"Length")
	<< "    e+ " << G4BestUnit(aCut->GetProductionCut("e+"),"Length")
	<< G4endl;
    log << " Energy thresholds : " ;
    log << " gamma " << G4BestUnit((*(table->GetEnergyCutsVector(0)))[aCouple->GetIndex()],"Energy")
	<< "    e- " << G4BestUnit((*(table->GetEnergyCutsVector(1)))[aCouple->GetIndex()],"Energy")
	<< "    e+ " << G4BestUnit((*(table->GetEnergyCutsVector(2)))[aCouple->GetIndex()],"Energy") << G4endl;
  } 
  G4cout << G4endl;

  log << "===================================================================" << G4endl;
  log<<"BremSplitting active ?: "<<ConfigData::BremSplitting::GetActivation()<<G4endl;
  if (ConfigData::BremSplitting::GetActivation()) {
    log<<"Factor                : "<<ConfigData::BremSplitting::GetFactor()<<G4endl;
  }
  log<<"# Secondaries         : "<<ConfigData::BremSplitting::GetNSecondaries()<<G4endl;
  log<<"Timing                : "<<*timer<<G4endl;

  log << "===================================================================" << G4endl;
  log<<"Number of events processed: "<<aRun->GetNumberOfEvent()<<G4endl;
}
