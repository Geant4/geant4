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
// $Id: PhotInRun.cc,v 1.6 2006/06/29 16:25:23 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

//#define debug

#include "PhotInRun.hh"

G4Allocator<PhotInRun> anPhotInRunAllocator;

PhotInRun::PhotInRun()
{
#ifdef debug
  G4cout<<"PhotInRun::PhotInRun Conctructor is called"<<G4endl;
#endif
  for(G4int i=0; i<PhotInDiNSections; i++)
  {
    totE[i] = 0.;
    totL[i] = 0.;
    nStep[i] = 0;
    nGamma[i] = 0;
    nElectron[i] = 0; 
    nPositron[i] = 0;
    eMinGamma[i] = DBL_MAX;
    eMinElectron[i] = DBL_MAX;
    eMinPositron[i] = DBL_MAX;
  }
}

PhotInRun::~PhotInRun() {}

void PhotInRun::RecordEvent(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
#ifdef debug
  G4cout<<"PhotInRun::RecordEvent: record event # "<<numberOfEvent<<G4endl;
#endif
  numberOfEvent++;
  for(G4int j=0; j<PhotInDiNSections; j++)
  {
    nGamma[j] += PhotInStackingAction::GetNGamma(j);
    nElectron[j] += PhotInStackingAction::GetNElectron(j);
    nPositron[j] += PhotInStackingAction::GetNPositron(j);
    if(eMinGamma[j]>PhotInStackingAction::GetEMinGamma(j))
    { eMinGamma[j] = PhotInStackingAction::GetEMinGamma(j); }
    if(eMinElectron[j]>PhotInStackingAction::GetEMinElectron(j))
    { eMinElectron[j] = PhotInStackingAction::GetEMinElectron(j); }
    if(eMinPositron[j]>PhotInStackingAction::GetEMinPositron(j))
    { eMinPositron[j] = PhotInStackingAction::GetEMinPositron(j); }
  }
  PhotInCalorHitsCollection* CHC = 0;
  for(G4int i=0; i<PhotInDiNSections; i++)
  {
    if (HCE) CHC = (PhotInCalorHitsCollection*)(HCE->GetHC(i));
    if (CHC)
    {
      G4int nHit = CHC->entries();
      for (G4int ii=0;ii<nHit;ii++)
      {
        totE[i] += (*CHC)[ii]->GetEDepos();
        totL[i] += (*CHC)[ii]->GetTrackL();
        nStep[i]+= (*CHC)[ii]->GetNSteps();
      }
    }
  }
#ifdef debug
  G4cout<<"PhotInRun::RecordEvent: End "<<G4endl;
#endif
}
