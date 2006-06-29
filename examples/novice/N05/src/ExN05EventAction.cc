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
// $Id: ExN05EventAction.cc,v 1.9 2006-06-29 17:53:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN05EventAction.hh"
#include "ExN05EventActionMessenger.hh"
#include "ExN05CalorimeterHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

ExN05EventAction::ExN05EventAction()
  :drawFlag(false)
{
  new ExN05EventActionMessenger(this);
}

ExN05EventAction::~ExN05EventAction()
{;}

void ExN05EventAction::BeginOfEventAction(const G4Event*)
{
  if(drawFlag)
    {
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
      if(pVVisManager)
	{
	  G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");
	}
    }
}

void ExN05EventAction::EndOfEventAction(const G4Event* evt )
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  calorimeterCollID    = SDman->GetCollectionID(colNam="CalCollection");
  hadCalorimeterCollID = SDman->GetCollectionID(colNam="HadCollection");
  
   G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  ExN05CalorimeterHitsCollection* CaloHC    = 0;
  ExN05CalorimeterHitsCollection* HadCaloHC = 0;
  if(HCE)
    {
      CaloHC    = (ExN05CalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
      HadCaloHC = (ExN05CalorimeterHitsCollection*)(HCE->GetHC(hadCalorimeterCollID));
    }
  
  if(CaloHC)
    {
      int n_hit = CaloHC->entries();
      G4cout << "     " << n_hit
	   << " hits are stored in EM ExN05CalorimeterHitsCollection." << G4endl;
      G4double totE = 0;
      for(int i=0;i<n_hit;i++)
	{ totE += (*CaloHC)[i]->GetEdep(); }
      G4cout << "     Total energy deposition in EM calorimeter crytals : "
	   << totE / GeV << " (GeV)" << G4endl;
    }

  if(HadCaloHC)
    {
      int n_hit = HadCaloHC->entries();
      G4cout << "     " << n_hit
	   << " hits are stored in HAD ExN05CalorimeterHitsCollection." << G4endl;
      G4double totE = 0;
      for(int i=0;i<n_hit;i++)
	{ totE += (*HadCaloHC)[i]->GetEdep(); }
      G4cout << "     Total energy deposition in HAD calorimeter towers : "
	   << totE / GeV << " (GeV)" << G4endl;
    }
  
  if(drawFlag)
    {
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
      if(pVVisManager)
	{
	  if(CaloHC)    CaloHC->DrawAllHits();
	  if(HadCaloHC) HadCaloHC->DrawAllHits();
	  G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
	}
    }
}



