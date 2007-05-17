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
// $Id: A01EventAction.cc,v 1.10 2007-05-17 09:55:14 duns Exp $
// --------------------------------------------------------------
//

#include "A01EventAction.hh"
#include "A01EventActionMessenger.hh"
#ifdef G4ANALYSIS_USE
#include "A01AnalysisManager.hh"
#endif // G4ANALYSIS_USE

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

#include "A01HodoscopeHit.hh"
#include "A01DriftChamberHit.hh"
#include "A01EmCalorimeterHit.hh"


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#include "A01HadCalorimeterHit.hh"
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

A01EventAction::A01EventAction()
{
  G4String colName;
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  HHC1ID = SDman->GetCollectionID(colName="hodoscope1/hodoscopeColl");
  HHC2ID = SDman->GetCollectionID(colName="hodoscope2/hodoscopeColl");
  DHC1ID = SDman->GetCollectionID(colName="chamber1/driftChamberColl");
  DHC2ID = SDman->GetCollectionID(colName="chamber2/driftChamberColl");
  ECHCID = SDman->GetCollectionID(colName="EMcalorimeter/EMcalorimeterColl");
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  HCHCID = SDman->GetCollectionID(colName="HadCalorimeter/HadCalorimeterColl");
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  verboseLevel = 1;
  messenger = new A01EventActionMessenger(this);

#ifdef G4ANALYSIS_USE
  plotter = 0;
  tuple = 0;
  dc1Hits = dc2Hits = 0;
  dc1XY = dc2XY = evstof = 0;

  // Do some analysis

  A01AnalysisManager* analysisManager = A01AnalysisManager::getInstance();
  IHistogramFactory* hFactory = analysisManager->getHistogramFactory();

  if (hFactory)
  {
    // Create some histograms
    dc1Hits = hFactory->createHistogram1D("Drift Chamber 1 # Hits",50,0,50);
    dc2Hits = hFactory->createHistogram1D("Drift Chamber 2 # Hits",50,0,50);

    // Create some clouds (Scatter Plots)
    dc1XY = hFactory->createCloud2D("Drift Chamber 1 X vs Y");
    dc2XY = hFactory->createCloud2D("Drift Chamber 2 X vs Y");
    evstof = hFactory->createCloud2D("EDep vs Time-of-flight");

    plotter = analysisManager->getPlotter();
    if (plotter)
    {
       plotter->createRegions(3,2);
       plotter->region(0)->plot(*dc1Hits);
       plotter->region(1)->plot(*dc2Hits);
       plotter->region(2)->plot(*dc1XY);
       plotter->region(3)->plot(*dc2XY);
       plotter->region(4)->plot(*evstof);
       plotter->show();
     }
  }

  // Create a Tuple

  ITupleFactory* tFactory = analysisManager->getTupleFactory();
  if (tFactory)
  {
     tuple = tFactory->create("MyTuple","MyTuple","int dc1Hits, dc2Hits, double ECEnergy, HCEnergy, time1, time2","");
  }
#endif // G4ANALYSIS_USE
}

A01EventAction::~A01EventAction()
{
#ifdef G4ANALYSIS_USE
  A01AnalysisManager::dispose();
#endif // G4ANALYSIS_USE
  delete messenger;
}

void A01EventAction::BeginOfEventAction(const G4Event*)
{
}

void A01EventAction::EndOfEventAction(const G4Event* evt)
{
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  A01HodoscopeHitsCollection* HHC1 = 0;
  A01HodoscopeHitsCollection* HHC2 = 0;
  A01DriftChamberHitsCollection* DHC1 = 0;
  A01DriftChamberHitsCollection* DHC2 = 0;
  A01EmCalorimeterHitsCollection* ECHC = 0;
  A01HadCalorimeterHitsCollection* HCHC = 0;
  if(HCE)
  {
    HHC1 = (A01HodoscopeHitsCollection*)(HCE->GetHC(HHC1ID));
    HHC2 = (A01HodoscopeHitsCollection*)(HCE->GetHC(HHC2ID));
    DHC1 = (A01DriftChamberHitsCollection*)(HCE->GetHC(DHC1ID));
    DHC2 = (A01DriftChamberHitsCollection*)(HCE->GetHC(DHC2ID));
    ECHC = (A01EmCalorimeterHitsCollection*)(HCE->GetHC(ECHCID));
    HCHC = (A01HadCalorimeterHitsCollection*)(HCE->GetHC(HCHCID));
  }

#ifdef G4ANALYSIS_USE
  // Fill some histograms

  if (DHC1 && dc1Hits)
  {
    int n_hit = DHC1->entries();
    dc1Hits->fill(n_hit);
    for(int i1=0;i1<n_hit;i1++)
    {
      A01DriftChamberHit* aHit = (*DHC1)[i1];
      G4ThreeVector localPos = aHit->GetLocalPos();
      if (dc1XY) dc1XY->fill(localPos.x(), localPos.y());
    }
  }
  if (DHC2 && dc2Hits)
  {
    int n_hit = DHC2->entries();
    dc2Hits->fill(n_hit);
    for(int i1=0;i1<n_hit;i1++)
    {
      A01DriftChamberHit* aHit = (*DHC2)[i1];
      G4ThreeVector localPos = aHit->GetLocalPos();
      if (dc2XY) dc2XY->fill(localPos.x(), localPos.y());
    }
  }

  // Fill the tuple

  if (tuple)
  {
	if (DHC1) tuple->fill(0,DHC1->entries());
	if (DHC2) tuple->fill(1,DHC2->entries());
	if(ECHC)
    {
      int iHit = 0;
      double totalE = 0.;
      for(int i1=0;i1<80;i1++)
      {
        A01EmCalorimeterHit* aHit = (*ECHC)[i1];
        double eDep = aHit->GetEdep();
        if(eDep>0.)
        {
          iHit++;
          totalE += eDep;
        }
      }
      tuple->fill(2,totalE);

	  if (HHC1 && HHC2 && HHC1->entries()==1 && HHC2->entries()==1)
	  {
	     double tof = (*HHC2)[0]->GetTime() - (*HHC1)[0]->GetTime();
		 if (evstof) evstof->fill(tof,totalE);
	  }
    }
    if(HCHC)
    {
      int iHit = 0;
      double totalE = 0.;
      for(int i1=0;i1<20;i1++)
      {
        A01HadCalorimeterHit* aHit = (*HCHC)[i1];
        double eDep = aHit->GetEdep();
        if(eDep>0.)
        {
          iHit++;
          totalE += eDep;
        }
      }
      tuple->fill(3,totalE);
    }
	if (HHC1 && HHC1->entries()==1) tuple->fill(4,(*HHC1)[0]->GetTime());
	if (HHC2 && HHC2->entries()==1) tuple->fill(5,(*HHC2)[0]->GetTime());
	tuple->addRow();
  }
  if (plotter) plotter->refresh();
#endif // G4ANALYSIS_USE


  // Diagnostics

  if (verboseLevel==0 || evt->GetEventID() % verboseLevel != 0) return;

  G4PrimaryParticle* primary = evt->GetPrimaryVertex(0)->GetPrimary(0);
  G4cout << G4endl
         << ">>> Event " << evt->GetEventID() << " >>> Simulation truth : "
         << primary->GetG4code()->GetParticleName()
         << " " << primary->GetMomentum() << G4endl;

  if(HHC1)
  {
    int n_hit = HHC1->entries();
    G4cout << "Hodoscope 1 has " << n_hit << " hits." << G4endl;
    for(int i1=0;i1<n_hit;i1++)
    {
      A01HodoscopeHit* aHit = (*HHC1)[i1];
      aHit->Print();
    }
  }
  if(HHC2)
  {
    int n_hit = HHC2->entries();
    G4cout << "Hodoscope 2 has " << n_hit << " hits." << G4endl;
    for(int i1=0;i1<n_hit;i1++)
    {
      A01HodoscopeHit* aHit = (*HHC2)[i1];
      aHit->Print();
    }
  }
  if(DHC1)
  {
    int n_hit = DHC1->entries();
    G4cout << "Drift Chamber 1 has " << n_hit << " hits." << G4endl;
    for(int i2=0;i2<5;i2++)
    {
      for(int i1=0;i1<n_hit;i1++)
      {
        A01DriftChamberHit* aHit = (*DHC1)[i1];
        if(aHit->GetLayerID()==i2) aHit->Print();
      }
    }
  }
  if(DHC2)
  {
    int n_hit = DHC2->entries();
    G4cout << "Drift Chamber 2 has " << n_hit << " hits." << G4endl;
    for(int i2=0;i2<5;i2++)
    {
      for(int i1=0;i1<n_hit;i1++)
      {
        A01DriftChamberHit* aHit = (*DHC2)[i1];
        if(aHit->GetLayerID()==i2) aHit->Print();
      }
    }
  }
  if(ECHC)
  {
    int iHit = 0;
    double totalE = 0.;
    for(int i1=0;i1<80;i1++)
    {
      A01EmCalorimeterHit* aHit = (*ECHC)[i1];
      double eDep = aHit->GetEdep();
      if(eDep>0.)
      {
        iHit++;
        totalE += eDep;
      }
    }
    G4cout << "EM Calorimeter has " << iHit << " hits. Total Edep is "
           << totalE/MeV << " (MeV)" << G4endl;
  }
  if(HCHC)
  {
    int iHit = 0;
    double totalE = 0.;
    for(int i1=0;i1<20;i1++)
    {
      A01HadCalorimeterHit* aHit = (*HCHC)[i1];
      double eDep = aHit->GetEdep();
      if(eDep>0.)
      {
        iHit++;
        totalE += eDep;
      }
    }
    G4cout << "Hadron Calorimeter has " << iHit << " hits. Total Edep is "
           << totalE/MeV << " (MeV)" << G4endl;
  }
}


