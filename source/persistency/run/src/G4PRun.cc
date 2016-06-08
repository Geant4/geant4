// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PRun.cc,v 1.1 1999/01/07 16:10:58 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#include "G4Run.hh"
#include "G4PRun.hh"

G4PRun::G4PRun()
// :runID(0),numberOfEvent(0),HCtable(NULL),DCtable(NULL)
:runID(0),numberOfEvent(0)
{;}

G4PRun::G4PRun(const G4Run* aRun)
{
  runID = aRun->GetRunID();
  numberOfEvent = aRun->GetNumberOfEvent();
}

G4PRun::~G4PRun()
{;}

