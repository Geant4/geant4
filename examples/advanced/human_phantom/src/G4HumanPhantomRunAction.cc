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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//

#include "G4HumanPhantomRunAction.hh"
#include "G4HumanPhantomAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Run.hh"

G4HumanPhantomRunAction::G4HumanPhantomRunAction(G4HumanPhantomAnalysisManager* analysis) :fAnalysisMan(analysis)
 { }

G4HumanPhantomRunAction::~G4HumanPhantomRunAction()
{}

void G4HumanPhantomRunAction::BeginOfRunAction(const G4Run* aRun)
{
 G4int run_number = aRun->GetRunID();
 G4cout << "### Run " << run_number << " start." << G4endl;
 fEnergyTotal["logicalHead"]=0.;
 fEnergyTotal["logicalTrunk"]=0.;
 fEnergyTotal["logicalLeftLeg"]=0.;
 fEnergyTotal["logicalRightLeg"]=0.;
 fEnergyTotal["logicalSkull"]=0.;
 fEnergyTotal["logicalLeftArmBone"]=0.;
 fEnergyTotal["logicalRightArmBone"]=0.;
 fEnergyTotal["logicalUpperSpine"]=0.;
 fEnergyTotal["logicalMiddleLowerSpine"]=0.;
 fEnergyTotal["logicalPelvis"]=0.;
 fEnergyTotal["logicalRibCage"]=0.;
 fEnergyTotal["logicalLeftClavicle"]=0.;
 fEnergyTotal["logicalRightClavicle"]=0.;
 fEnergyTotal["logicalLeftLegBone"]=0.;
 fEnergyTotal["logicalRightLegBone"]=0.;
 fEnergyTotal["logicalLeftScapula"]=0.; 
 fEnergyTotal["logicalRightScapula"]=0.;
 fEnergyTotal["logicalHeart"]=0.;
 fEnergyTotal["logicalThyroid"]=0.;
 fEnergyTotal["logicalThymus"]=0.;
 fEnergyTotal["logicalMaleGenitalia"]=0.;
 fEnergyTotal["logicalBrain"]=0.;
 fEnergyTotal["logicalStomach"]=0.;
 fEnergyTotal["logicalUpperLargeIntestine"]=0.;
 fEnergyTotal["logicalLowerLargeIntestine"]=0.;
 fEnergyTotal["logicalSmallIntestine"]=0;
 fEnergyTotal["logicalSpleen"]=0.;
 fEnergyTotal["logicalPancreas"]=0.;
 fEnergyTotal["logicalLeftKidney"]=0.;
 fEnergyTotal["logicalRightKidney"]=0.;
 fEnergyTotal["logicalUrinaryBladder"]=0.;
 fEnergyTotal["logicalUterus"]=0.;
 fEnergyTotal["logicalLeftLung"]=0.;
 fEnergyTotal["logicalRightLung"]=0.;
 fEnergyTotal["logicalLeftOvary"]=0.;
 fEnergyTotal["logicalRightOvary"]=0.;
 fEnergyTotal["logicalLeftTeste"]=0;
 fEnergyTotal["logicalRightTeste"]=0;
 fEnergyTotal["logicalLeftBreast"]=0.;
 fEnergyTotal["logicalRightBreast"]=0.; 
 fEnergyTotal["logicalLeftAdrenal"]=0.; 
 fEnergyTotal["logicalRightAdrenal"]=0.;

 // Create ROOT file, histograms and ntuple
 fAnalysisMan -> book();
}

void G4HumanPhantomRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Number of events = " << aRun->GetNumberOfEvent() << G4endl;
  totalRunEnergyDeposit();

// Close the output ROOT file with the results
   fAnalysisMan -> save(); 
}

void G4HumanPhantomRunAction::Fill(G4String bName, 
				      G4double energyDeposit)
{
 fEnergyTotal[bName] += energyDeposit;
}

void G4HumanPhantomRunAction::totalRunEnergyDeposit() 
{
 std::map<std::string,G4double>::iterator i = fEnergyTotal.begin();
  std::map<std::string,G4double>::iterator end = fEnergyTotal.end();

  G4double totalEnergyDepositInPhantom =0.;
  G4int k=0;
  while(i!=end)
    {
      G4String bodypart = i->first;
      G4double energyDep = i->second;
      //  if(energyDep != 0.)
      //	{
     
      G4cout << "Energy Total in Run:" << bodypart << ", ID: "  << k 
             << ", Energy Deposition (MeV): "  
	     << energyDep/MeV
	     << G4endl;

     // Fill Ntuple
     fAnalysisMan -> FillNtupleWithEnergyDeposition(k, energyDep/MeV);

     i++;
     k++;

     totalEnergyDepositInPhantom += energyDep;
    }
  
 G4cout << "Total Energy deposit in the body is: " 
	<< totalEnergyDepositInPhantom/MeV << " MeV" <<G4endl; 
}
