// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05RunAction.hh,v 1.2 1999-04-16 12:04:59 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN05RunAction_h
#define ExN05RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"

class G4Run;

class ExN05RunAction : public G4UserRunAction
{
  public:
    ExN05RunAction();
    ~ExN05RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

  private:
    G4Timer* timer;
    G4int runIDcounter;
};

#endif

