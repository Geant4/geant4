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
// $Id: G4PRun.cc,v 1.7 2001/07/11 10:02:27 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "G4Run.hh"
#include "G4PRun.hh"
#include "G4Event.hh"

G4PRun::G4PRun()
// :runID(0),numberOfEvent(0),HCtable(0),DCtable(0)
:runID(0),numberOfEvent(0)
{;}

G4PRun::G4PRun(const G4Run* aRun)
{
  runID = aRun->GetRunID();
  numberOfEvent = aRun->GetNumberOfEvent();
}

G4PRun::~G4PRun()
{;}

G4Run* G4PRun::MakeTransientObject()
{
  G4Run* aRun = new G4Run();

  aRun->SetRunID(runID);

  G4Event* dummyEvt = new G4Event;

  for(G4int i = 0; i < numberOfEvent; i++)
  {
    aRun->RecordEvent(dummyEvt);
  }

  delete dummyEvt;

  return aRun;
}

