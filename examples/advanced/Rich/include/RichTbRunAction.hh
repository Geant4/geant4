// Rich advanced example for Geant4
// RichTbRunAction.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbRunAction_h
#define RichTbRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class RichTbRunAction : public G4UserRunAction
{
  public:

  RichTbRunAction();

  virtual ~RichTbRunAction();

  public:

  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);

  private:

    G4Timer* timer;
};

#endif 


