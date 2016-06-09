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
// $Id: ExN07Run.cc,v 1.3 2003/04/09 23:20:59 asaim Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#include "ExN07Run.hh"
#include "ExN07CalorHit.hh"
#include "ExN07StackingAction.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"

G4Allocator<ExN07Run> anExN07RunAllocator;

ExN07Run::ExN07Run()
{
  for(size_t i=0;i<6;i++)
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

ExN07Run::~ExN07Run()
{;}

void ExN07Run::RecordEvent(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  for(int j=0;j<6;j++)
  {
    nGamma[j] += ExN07StackingAction::GetNGamma(j);
    nElectron[j] += ExN07StackingAction::GetNElectron(j);
    nPositron[j] += ExN07StackingAction::GetNPositron(j);
    if(eMinGamma[j]>ExN07StackingAction::GetEMinGamma(j))
    { eMinGamma[j] = ExN07StackingAction::GetEMinGamma(j); }
    if(eMinElectron[j]>ExN07StackingAction::GetEMinElectron(j))
    { eMinElectron[j] = ExN07StackingAction::GetEMinElectron(j); }
    if(eMinPositron[j]>ExN07StackingAction::GetEMinPositron(j))
    { eMinPositron[j] = ExN07StackingAction::GetEMinPositron(j); }
  }
  ExN07CalorHitsCollection* CHC = 0;
  for(size_t i=0;i<6;i++)
  {
    if (HCE) CHC = (ExN07CalorHitsCollection*)(HCE->GetHC(i));
    if (CHC)
    {
      G4int nHit = CHC->entries();
      for (G4int ii=0;ii<nHit;ii++)
      {
        totE[i] += (*CHC)[ii]->GetEdep();
        totL[i] += (*CHC)[ii]->GetTrak();
        nStep[i] += (*CHC)[ii]->GetNStep();
      }
    }
  }
}

  
