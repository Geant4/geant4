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
/// \file analysis/A01/src/A01EventAction.cc
/// \brief Implementation of the A01EventAction class
//
// $Id$
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
#include "G4SystemOfUnits.hh"

#include "A01HodoscopeHit.hh"
#include "A01DriftChamberHit.hh"
#include "A01EmCalorimeterHit.hh"

#include "A01HadCalorimeterHit.hh"

A01EventAction::A01EventAction()
{
  G4String colName;
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  fHHC1ID = SDman->GetCollectionID(colName="hodoscope1/hodoscopeColl");
  fHHC2ID = SDman->GetCollectionID(colName="hodoscope2/hodoscopeColl");
  fDHC1ID = SDman->GetCollectionID(colName="chamber1/driftChamberColl");
  fDHC2ID = SDman->GetCollectionID(colName="chamber2/driftChamberColl");
  fECHCID = SDman->GetCollectionID(colName="EMcalorimeter/EMcalorimeterColl");
  fHCHCID = SDman->GetCollectionID(colName="HadCalorimeter/HadCalorimeterColl");
  fVerboseLevel = 1;
  fMessenger = new A01EventActionMessenger(this);

#ifdef G4ANALYSIS_USE
  fPlotter = 0;
  fTuple = 0;
  fDc1Hits = fDc2Hits = 0;
  fDc1XY = fDc2XY = fEvstof = 0;

  // Do some analysis

  A01AnalysisManager* analysisManager = A01AnalysisManager::getInstance();
  IHistogramFactory* hFactory = analysisManager->getHistogramFactory();

  if (hFactory)
  {
    // Create some histograms
    fDc1Hits = hFactory->createHistogram1D("Drift Chamber 1 # Hits",50,0,50);
    fDc2Hits = hFactory->createHistogram1D("Drift Chamber 2 # Hits",50,0,50);

    // Create some clouds (Scatter Plots)
    fDc1XY = hFactory->createCloud2D("Drift Chamber 1 X vs Y");
    fDc2XY = hFactory->createCloud2D("Drift Chamber 2 X vs Y");
    fEvstof = hFactory->createCloud2D("EDep vs Time-of-flight");

    fPlotter = analysisManager->getPlotter();
    if (fPlotter)
    {
       fPlotter->createRegions(3,2);
       fPlotter->region(0)->plot(*fDc1Hits);
       fPlotter->region(1)->plot(*fDc2Hits);
       fPlotter->region(2)->plot(*fDc1XY);
       fPlotter->region(3)->plot(*fDc2XY);
       fPlotter->region(4)->plot(*fEvstof);
       fPlotter->show();
     }
  }

  // Create a Tuple

  ITupleFactory* tFactory = analysisManager->getTupleFactory();
  if (tFactory)
  {
     fTuple = tFactory->create("MyTuple","MyTuple","int fDc1Hits, fDc2Hits, double ECEnergy, HCEnergy, time1, time2","");
  }
#endif // G4ANALYSIS_USE
}

A01EventAction::~A01EventAction()
{
#ifdef G4ANALYSIS_USE
  A01AnalysisManager::dispose();
#endif // G4ANALYSIS_USE
  delete fMessenger;
}

void A01EventAction::BeginOfEventAction(const G4Event*)
{
}

void A01EventAction::EndOfEventAction(const G4Event* evt)
{
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  A01HodoscopeHitsCollection* fHHC1 = 0;
  A01HodoscopeHitsCollection* fHHC2 = 0;
  A01DriftChamberHitsCollection* fDHC1 = 0;
  A01DriftChamberHitsCollection* fDHC2 = 0;
  A01EmCalorimeterHitsCollection* ECHC = 0;
  A01HadCalorimeterHitsCollection* HCHC = 0;
  if(HCE)
  {
    fHHC1 = (A01HodoscopeHitsCollection*)(HCE->GetHC(fHHC1ID));
    fHHC2 = (A01HodoscopeHitsCollection*)(HCE->GetHC(fHHC2ID));
    fDHC1 = (A01DriftChamberHitsCollection*)(HCE->GetHC(fDHC1ID));
    fDHC2 = (A01DriftChamberHitsCollection*)(HCE->GetHC(fDHC2ID));
    ECHC = (A01EmCalorimeterHitsCollection*)(HCE->GetHC(fECHCID));
    HCHC = (A01HadCalorimeterHitsCollection*)(HCE->GetHC(fHCHCID));
  }

#ifdef G4ANALYSIS_USE
  // Fill some histograms

  if (fDHC1 && fDc1Hits)
  {
    int n_hit = fDHC1->entries();
    fDc1Hits->fill(n_hit);
    for(int i1=0;i1<n_hit;i1++)
    {
      A01DriftChamberHit* aHit = (*fDHC1)[i1];
      G4ThreeVector localPos = aHit->GetLocalPos();
      if (fDc1XY) fDc1XY->fill(localPos.y(), localPos.x());
    }
  }
  if (fDHC2 && fDc2Hits)
  {
    int n_hit = fDHC2->entries();
    fDc2Hits->fill(n_hit);
    for(int i1=0;i1<n_hit;i1++)
    {
      A01DriftChamberHit* aHit = (*fDHC2)[i1];
      G4ThreeVector localPos = aHit->GetLocalPos();
      if (fDc2XY) fDc2XY->fill(localPos.y(), localPos.x());
    }
  }

  // Fill the tuple

  if (fTuple)
  {
    if (fDHC1) fTuple->fill(0,fDHC1->entries());
    if (fDHC2) fTuple->fill(1,fDHC2->entries());
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
      fTuple->fill(2,totalE);

          if (fHHC1 && fHHC2 && fHHC1->entries()==1 && fHHC2->entries()==1)
          {
             double tof = (*fHHC2)[0]->GetTime() - (*fHHC1)[0]->GetTime();
                 if (fEvstof) fEvstof->fill(tof,totalE);
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
      fTuple->fill(3,totalE);
    }
    if (fHHC1 && fHHC1->entries()==1) fTuple->fill(4,(*fHHC1)[0]->GetTime());
    if (fHHC2 && fHHC2->entries()==1) fTuple->fill(5,(*fHHC2)[0]->GetTime());
        fTuple->addRow();
  }
  if (fPlotter) fPlotter->refresh();
#endif // G4ANALYSIS_USE


  // Diagnostics

  if (fVerboseLevel==0 || evt->GetEventID() % fVerboseLevel != 0) return;

  G4PrimaryParticle* primary = evt->GetPrimaryVertex(0)->GetPrimary(0);
  G4cout << G4endl
         << ">>> Event " << evt->GetEventID() << " >>> Simulation truth : "
         << primary->GetG4code()->GetParticleName()
         << " " << primary->GetMomentum() << G4endl;

  if(fHHC1)
  {
    int n_hit = fHHC1->entries();
    G4cout << "Hodoscope 1 has " << n_hit << " hits." << G4endl;
    for(int i1=0;i1<n_hit;i1++)
    {
      A01HodoscopeHit* aHit = (*fHHC1)[i1];
      aHit->Print();
    }
  }
  if(fHHC2)
  {
    int n_hit = fHHC2->entries();
    G4cout << "Hodoscope 2 has " << n_hit << " hits." << G4endl;
    for(int i1=0;i1<n_hit;i1++)
    {
      A01HodoscopeHit* aHit = (*fHHC2)[i1];
      aHit->Print();
    }
  }
  if(fDHC1)
  {
    int n_hit = fDHC1->entries();
    G4cout << "Drift Chamber 1 has " << n_hit << " hits." << G4endl;
    for(int i2=0;i2<5;i2++)
    {
      for(int i1=0;i1<n_hit;i1++)
      {
        A01DriftChamberHit* aHit = (*fDHC1)[i1];
        if(aHit->GetLayerID()==i2) aHit->Print();
      }
    }
  }
  if(fDHC2)
  {
    int n_hit = fDHC2->entries();
    G4cout << "Drift Chamber 2 has " << n_hit << " hits." << G4endl;
    for(int i2=0;i2<5;i2++)
    {
      for(int i1=0;i1<n_hit;i1++)
      {
        A01DriftChamberHit* aHit = (*fDHC2)[i1];
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


