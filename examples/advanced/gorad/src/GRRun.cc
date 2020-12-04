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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRRun.cc
//   Gorad Run class that handles filling histograms and profile plots
//   with scores accumulated by scoeres for each event.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRRun.hh"
#include "GRRunAction.hh"
#include "GRAnalysis.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

GRRun::GRRun(GRRunAction* ra) : G4Run(),pRA(ra)
{
  ;
}

GRRun::~GRRun()
{
  ;
}

void GRRun::RecordEvent(const G4Event* anEvent)
{
  
  numberOfEvent++;  // This is an original line.

  G4HCofThisEvent* pHCE = anEvent->GetHCofThisEvent();

  auto analysisManager = G4AnalysisManager::Instance();
  auto map = pRA->IDMap;
  for(auto itr : map)
  {
    if(itr.second->pplotter!=nullptr) continue; // directly plotted by the PrimitivePlotter

    auto cID = itr.second->collID;
    auto hID = itr.second->histID;
    auto hTyp = itr.second->histType;

    if(hTyp==1) // 1D histogram
    {
      if(cID>=0) // scorer
      {
        if(!pHCE) continue;
        auto score = (G4THitsMap<G4double>*)(pHCE->GetHC(cID));
        G4double val = 0.;
        for(auto hItr : *score)
        {
          if(itr.second->idx==-1 || itr.second->idx==hItr.first)
          { val += *(hItr.second); }
        }
        analysisManager->FillH1(hID,val);
      }
      else // primary particle
      {
        auto pv = anEvent->GetPrimaryVertex();
        while(pv)
        {
          auto pp = pv->GetPrimary();
          while(pp)
          {
            auto primE = pp->GetKineticEnergy();
            G4double weight = 1.0;
            if(itr.second->biasf) weight = pp->GetWeight();
            analysisManager->FillH1(hID,primE,weight);
            pp = pp->GetNext();
          }
          pv = pv->GetNext();
        }
      }  
    }

    else if(hTyp==2) // 1D profile plot
    {
      if(!pHCE) continue;
      auto score = (G4THitsMap<G4double>*)(pHCE->GetHC(cID));
      for(auto hItr : *score)
      { analysisManager->FillP1(hID,G4double(hItr.first),*(hItr.second)); }
    }
    
  }

  auto ntmap = pRA->NTMap;
  if(ntmap.size()>0) 
  {
    for(auto ntitr : ntmap)
    {
      auto colID = ntitr.first;
      auto cID = ntitr.second->collID;
      G4double val = 0.;
      if(cID>=0)
      {
        if(!pHCE) continue;
        auto score = (G4THitsMap<G4double>*)(pHCE->GetHC(cID));
        for(auto hItr : *score)
        {
          if(ntitr.second->idx==-1 || ntitr.second->idx==hItr.first)
          { val += *(hItr.second); }
        }
      }
      else
      {
        auto pv = anEvent->GetPrimaryVertex();
        auto pp = pv->GetPrimary();
        val = pp->GetKineticEnergy();
      }
      analysisManager->FillNtupleDColumn(colID,val*(ntitr.second->fuct));
    }
    analysisManager->AddNtupleRow();
  }

}

void GRRun::Merge(const G4Run * aRun) 
{
  G4Run::Merge(aRun);
}


