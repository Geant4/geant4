// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVARunAction.hh,v 1.3 2001-02-07 17:31:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVARunAction_h
#define TstVARunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Run;

class TstVARunAction : public G4UserRunAction
{
  public:
    TstVARunAction();
    virtual ~TstVARunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

