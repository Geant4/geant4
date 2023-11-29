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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//

#include <stdlib.h>
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "AnalysisMessenger.hh"

AnalysisManager::AnalysisManager(AnalysisMessenger* analysisMessenger) 
{
  factoryOn = false;

// Initialization
// histograms
  //for (G4int k=0; k<MaxHisto; k++) fHistId[k] = 0;

// Initialization ntuple
  for (G4int k=0; k<MaxNtCol; k++) fNtColId[k] = 0;

  //h10 = 0;
  //h20 = 0;
  
  messenger = analysisMessenger;
}

AnalysisManager::~AnalysisManager() 
{
}

void AnalysisManager::book(G4bool addExtraNt) 
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  
  manager->SetVerboseLevel(2);
 
  extraNt = addExtraNt;
  
  usingRoot = messenger -> IsRootOutput();
  
  // Create an output file
  G4String fileName;
  if( usingRoot == true ) fileName = "radioprotection.root";
  else fileName = "radioprotection.csv";

  // Create directories (not supported by csv)
  if( usingRoot == true ) manager->SetNtupleDirectoryName("radioprotection_ntuple");
  

  G4bool fileOpen = manager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fileName << G4endl;
    return;
  }

  manager->SetFirstNtupleId(1);

  //Create Primary Energy Ntuple
  manager -> CreateNtuple("101", "Primary Energy");
  fNtColId[0] = manager -> CreateNtupleDColumn("Ek");
  manager -> FinishNtuple();

  //Create Energy Deposition and Path Length within SV Ntuple
  manager -> CreateNtuple("102", "Edep");
  fNtColId[1] = manager -> CreateNtupleDColumn("edep");
  fNtColId[2] = manager -> CreateNtupleDColumn("len");
  if( extraNt == true ) fNtColId[3] = manager -> CreateNtupleIColumn("eventID");
  manager -> FinishNtuple();

  //creating a ntuple, containing the information about secondary particles
  manager -> CreateNtuple("103", "secondary");
  fNtColId[4] = manager -> CreateNtupleDColumn("AA");
  fNtColId[5] = manager -> CreateNtupleDColumn("ZZ");
  fNtColId[6] = manager -> CreateNtupleDColumn("KE");
  manager -> FinishNtuple();
  
  // if the telescope is in use, add an exta ntuple for its information
  if( extraNt == true)
  {
      manager -> CreateNtuple("104", "secondStageE");
      fNtColId[7] = manager -> CreateNtupleDColumn("edep");
      fNtColId[8] = manager -> CreateNtupleIColumn("eventID");
      manager -> FinishNtuple();
  }
  
  factoryOn = true;    
}


void AnalysisManager::SetPrimaryEnergy(G4double energy)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(1, fNtColId[0], energy);
  manager -> AddNtupleRow(1); 
}

void AnalysisManager::StoreEnergyDeposition(G4double edep, G4double len, G4int eid)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(2, fNtColId[1], edep);
  manager -> FillNtupleDColumn(2, fNtColId[2], len);
  if( extraNt == true ) manager -> FillNtupleIColumn(2, fNtColId[3], eid);
    // the event ID is only useful when the 4th ntuple is enabled
  manager -> AddNtupleRow(2); 
}

void AnalysisManager::FillSecondaries(G4int AA, G4double charge, G4double energy)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(3, fNtColId[4], AA);
  manager -> FillNtupleDColumn(3, fNtColId[5], charge);
  manager -> FillNtupleDColumn(3, fNtColId[6], energy);
  manager -> AddNtupleRow(3);  
}

void AnalysisManager::StoreSecondStageEnergyDeposition(G4double edep, G4int eid)
{
  if( extraNt == true )
  {
    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager -> FillNtupleDColumn(4, fNtColId[7], edep);
    manager -> FillNtupleIColumn(4, fNtColId[8], eid);
    manager -> AddNtupleRow(4);
  }
}

void AnalysisManager::finish() 
{   
 if (factoryOn) 
   {
    G4AnalysisManager* manager = G4AnalysisManager::Instance();    
    manager -> Write();
    manager -> CloseFile();  
    factoryOn = false;
   }
}













