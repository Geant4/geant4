// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02RunAction.hh,v 1.3 1999/11/29 18:23:32 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#ifndef PersEx02RunAction_h
#define PersEx02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class PersEx02RunAction : public G4UserRunAction
{
  public:
    PersEx02RunAction();
    virtual ~PersEx02RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
};

#endif

