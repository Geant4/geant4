// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserRunAction.hh,v 1.1 1999-01-07 16:14:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UserRunAction_h
#define G4UserRunAction_h 1

class G4Run;

class G4UserRunAction
{
  public:
    G4UserRunAction();
    virtual ~G4UserRunAction();

  public:
    virtual void BeginOfRunAction(G4Run* aRun);
    virtual void EndOfRunAction(G4Run* aRun);
};

#endif


