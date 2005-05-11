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
// $Id: PhotInRun.cc,v 1.1 2005-05-11 10:37:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "PhotInRun.hh"
#include "PhotInCalorHit.hh"
#include "PhotInStackingAction.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"

G4Allocator<PhotInRun> anPhotInRunAllocator;

PhotInRun::PhotInRun()
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

PhotInRun::~PhotInRun()
{;}

void PhotInRun::RecordEvent(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  for(int j=0;j<6;j++)
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
  for(size_t i=0;i<6;i++)
  {
    if (HCE) CHC = (PhotInCalorHitsCollection*)(HCE->GetHC(i));
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

  
