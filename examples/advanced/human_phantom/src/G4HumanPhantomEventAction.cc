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
hitCollectionID(-1)
{ 
 
}
 
G4HumanPhantomEventAction::~G4HumanPhantomEventAction()
{
}

void G4HumanPhantomEventAction::BeginOfEventAction(const G4Event*)
{
 energyTotal["logicalHead"]=0.;
 energyTotal["logicalTrunk"]=0.;
 energyTotal["logicalLeftLeg"]=0.;
 energyTotal["logicalRightLeg"]=0.;
 energyTotal["logicalSkull"]=0.;
 energyTotal["logicalLeftArmBone"]=0.;
 energyTotal["logicalRightArmBone"]=0.;
 energyTotal["logicalUpperSpine"]=0.;
 energyTotal["logicalMiddleLowerSpine"]=0.;
 energyTotal["logicalPelvis"]=0.;
 energyTotal["logicalRibCage"]=0.;
 energyTotal["logicalLeftClavicle"]=0.;
 energyTotal["logicalRightClavicle"]=0.;
 energyTotal["logicalLeftLegBone"]=0.;
 energyTotal["logicalRightLegBone"]=0.;
 energyTotal["logicalLeftScapula"]=0.; 
 energyTotal["logicalRightScapula"]=0.;
 energyTotal["logicalHeart"]=0.;
 energyTotal["logicalThyroid"]=0.;
 energyTotal["logicalThymus"]=0.;
 energyTotal["logicalMaleGenitalia"]=0.;
 energyTotal["logicalBrain"]=0.;
 energyTotal["logicalStomach"]=0.;
 energyTotal["logicalUpperLargeIntestine"]=0.;
 energyTotal["logicalLowerLargeIntestine"]=0.;
 energyTotal["logicalSmallIntestine"]=0;
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
 energyTotal["logicalLeftTeste"]=0;
 energyTotal["logicalRightTeste"]=0;
 energyTotal["logicalLeftBreast"]=0.;
 energyTotal["logicalRightBreast"]=0.; 
 energyTotal["logicalLeftAdrenal"]=0.; 
 energyTotal["logicalRightAdrenal"]=0.;

 G4SDManager * SDman = G4SDManager::GetSDMpointer();  

  if (hitCollectionID==-1) {
    hitCollectionID = SDman->GetCollectionID("HumanPhantomCollection");
  }
}
 
void G4HumanPhantomEventAction::EndOfEventAction(const G4Event* evt)
{  

 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 
 G4HumanPhantomHitsCollection* HC = 0;

 if (HCE)
     HC = (G4HumanPhantomHitsCollection*)(HCE->GetHC(hitCollectionID));

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
 energyTotal[bName] += energyDeposit;
}

void G4HumanPhantomEventAction::totalEventEnergyDeposit() 
{

 G4RunManager* runManager = G4RunManager::GetRunManager();
 G4HumanPhantomRunAction* pointerRun = (G4HumanPhantomRunAction*)(runManager->GetUserRunAction());

 std::map<std::string,G4double>::iterator i = energyTotal.begin();
  std::map<std::string,G4double>::iterator end = energyTotal.end();

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
