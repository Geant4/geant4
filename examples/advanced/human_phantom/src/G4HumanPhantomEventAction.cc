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
#include "G4HumanPhantomEventAction.hh"
#include "G4HumanPhantomHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4HumanPhantomRunAction.hh"
#include "G4RunManager.hh"

G4HumanPhantomEventAction::G4HumanPhantomEventAction():
fHitCollectionID(-1)
{ 
 
}

void G4HumanPhantomEventAction::BeginOfEventAction(const G4Event*)
{
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

 G4SDManager * SDman = G4SDManager::GetSDMpointer();  

  if (fHitCollectionID==-1) {
    fHitCollectionID = SDman->GetCollectionID("HumanPhantomCollection");
  }
}
 
void G4HumanPhantomEventAction::EndOfEventAction(const G4Event* evt)
{  

 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 
 G4HumanPhantomHitsCollection* HC = nullptr;

 if (HCE)
     HC = (G4HumanPhantomHitsCollection*)(HCE->GetHC(fHitCollectionID));

 if (HC)
	{
	  G4int hitNumber = HC->entries();
          G4double edep =0;
          G4String bodyPart;
          for (G4int i=0;i<hitNumber;i++) 
	    {
	      edep = (*HC)[i]->GetEdep();
	      bodyPart = (*HC)[i]->GetBodyPartID();
              Fill(bodyPart, edep);
	      //	      if(edep !=0.) G4cout << bodyPart <<": "<< edep/MeV << G4endl; 
	    }
	}

 totalEventEnergyDeposit();
}

void G4HumanPhantomEventAction:: Fill(G4String bName, 
				      G4double energyDeposit)

{
 fEnergyTotal[bName] += energyDeposit;
}

void G4HumanPhantomEventAction::totalEventEnergyDeposit() 
{

 G4RunManager* runManager = G4RunManager::GetRunManager();
 auto* pointerRun = (G4HumanPhantomRunAction*)(runManager->GetUserRunAction());

 std::map<std::string,G4double>::iterator i = fEnergyTotal.begin();
  std::map<std::string,G4double>::iterator end = fEnergyTotal.end();

  while(i!=end)
    {

      G4String bodypart = i->first;
      G4double energyDep = i->second;
      
      if(energyDep != 0.)
	{
	  pointerRun->Fill(bodypart, energyDep);
	}
      i++;
    }
  
}
