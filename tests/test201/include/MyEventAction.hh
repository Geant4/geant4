// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyEventAction.hh,v 1.3 1999-12-15 14:54:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyEventAction_h
#define MyEventAction_h 1

#include "G4UserEventAction.hh"

class G4Event;

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction();
    ~MyEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
};

#endif

    
