// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08EventAction.hh,v 1.1 1999-01-08 16:35:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08EventAction_h
#define T08EventAction_h 1

#include "G4UserEventAction.hh"

class G4Event;

class T08EventAction : public G4UserEventAction
{
  public:
    T08EventAction();
    ~T08EventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();
    
  private:

};

#endif

    
