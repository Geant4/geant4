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
// $Id: ExN07Run.cc,v 1.1 2003-03-10 01:43:37 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN07Run.hh"
#include "ExN07CalorHit.hh"

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
  }
}

ExN07Run::~ExN07Run()
{;}

void ExN07Run::RecordEvent(G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
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

  
