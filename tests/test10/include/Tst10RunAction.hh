// $Id: Tst10RunAction.hh,v 1.1 1999-01-08 16:35:33 gunter Exp $
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
    ~Tst10RunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

