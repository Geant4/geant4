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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#include "G4HumanPhantomRunAction.hh"
#include "G4HumanPhantomAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Run.hh"

G4HumanPhantomRunAction::G4HumanPhantomRunAction()
{
 // Create analysis manager
 G4AnalysisManager::Instance();
}

G4HumanPhantomRunAction::~G4HumanPhantomRunAction()
{
delete G4AnalysisManager::Instance();
}

void G4HumanPhantomRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int run_number = aRun->GetRunID();
  G4cout << "### Run " << run_number << " start." << G4endl;

 energyTotal["logicalHead"]=0.;
 energyTotal["logicalTrunk"]=0.;
 energyTotal["logicalLeftLeg"]=0.;
 energyTotal["logicalRightLeg"]=0.;
 energyTotal["logicalBrain"]=0.;
 energyTotal["logicalLeftArmBone"]=0.;
 energyTotal["logicalRightArmBone"]=0.;
 energyTotal["logicalSkull"]=0.;
 energyTotal["logicalUpperSpine"]=0.;
 energyTotal["logicalMiddleLowerSpine"]=0.;
 energyTotal["logicalPelvis"]=0.;
 energyTotal["logicalStomach"]=0.;
 energyTotal["logicalUpperLargeIntestine"]=0.;
 energyTotal["logicalLowerLargeIntestine"]=0.;
 energyTotal["logicalRibCage"]=0.;
 energyTotal["logicalSpleen"]=0.;
 energyTotal["logicalPancreas"]=0.;
 energyTotal["logicalLeftKidney"]=0.;
 energyTotal["logicalRightKidney"]=0.;
 energyTotal["logicalUrinaryBladder"]=0.;
 energyTotal["logicalUterus"]=0.;
 energyTotal["logicalLeftLung"]=0.;
 energyTotal["logicalRightLung"]=0.;
 energyTotal["logicalLeftOvary"]=0.;
 energyTotal["logicalRightOvary"]=0.;
 energyTotal["logicalLeftLegBone"]=0.;
 energyTotal["logicalRightLegBone"]=0.;
 energyTotal["logicalLeftBreast"]=0.;
 energyTotal["logicalRightBreast"]=0.; 
 energyTotal["logicalLeftScapula"]=0.; 
 energyTotal["logicalRightScapula"]=0.; 
 energyTotal["logicalLeftAdrenal"]=0.; 
 energyTotal["logicalRightAdrenal"]=0.; 

// Create ROOT output file with ntuple containing the
// energy deposition in the organs.

// Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Open a ROOT output file
  //
  G4String fileName = "g4humanphantom.root";
  analysisManager -> OpenFile(fileName);
 
  // Create the ntuple
  analysisManager -> CreateNtuple("ntuple1", "Edep (MeV) in the organs");
  analysisManager -> CreateNtupleDColumn("ID"); // Integer identifying the organ
  analysisManager -> CreateNtupleDColumn("Edep"); // Energy deposition expressed in MeV
  analysisManager -> FinishNtuple();
}

void G4HumanPhantomRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Number of events = " << aRun->GetNumberOfEvent() << G4endl;
  totalRunEnergyDeposit();

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Save the ntuple and histograms and close the ROOT file
  analysisManager -> Write();
  analysisManager -> CloseFile();
}


void G4HumanPhantomRunAction::Fill(G4String bName, 
				      G4double energyDeposit)

{
 energyTotal[bName] += energyDeposit;
}

void G4HumanPhantomRunAction::totalRunEnergyDeposit() 
{
 std::map<std::string,G4double>::iterator i = energyTotal.begin();
  std::map<std::string,G4double>::iterator end = energyTotal.end();

  G4double totalEnergyDepositInPhantom =0.;
  G4int k=0;
  while(i!=end)
    {
      G4String bodypart = i->first;
      G4double energyDep = i->second;
      //  if(energyDep != 0.)
      //	{
     
      G4cout << "Energy Total in Run" <<bodypart << " = "  
	     << G4BestUnit(energyDep,"Energy") 
	     << G4endl;

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Fill ntuple
  analysisManager -> FillNtupleDColumn(0, k);
  analysisManager -> FillNtupleDColumn(1, energyDep/MeV);
  analysisManager -> AddNtupleRow();  
  //	}
      i++;
      k++;
      totalEnergyDepositInPhantom += energyDep;
    }
  
 
 G4cout << "Total Energy deposit in the body is: " 
	<< totalEnergyDepositInPhantom/MeV << " MeV" <<G4endl; 
}



