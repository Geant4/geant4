// $Id: Tst10RunAction.hh,v 1.2 1999-04-17 08:01:49 kurasige Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserRunList
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Tst10RunAction_h
#define Tst10RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst10RunAction : public G4UserRunAction
{
  public:
    Tst10RunAction();
    virtual ~Tst10RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

