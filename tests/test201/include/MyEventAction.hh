// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyEventAction.hh,v 1.1 1999-01-08 16:35:53 gunter Exp $
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
    void BeginOfEventAction();
    void EndOfEventAction();
};

#endif

    
