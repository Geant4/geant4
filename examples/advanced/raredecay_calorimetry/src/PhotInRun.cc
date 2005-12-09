//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PhotInRun.cc,v 1.5 2005-12-09 16:44:21 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
