// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02EventAction.hh,v 1.1 1999-01-07 16:05:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef ExN02EventAction_h
#define ExN02EventAction_h 1

#include "G4UserEventAction.hh"

class G4Event;

class ExN02EventAction : public G4UserEventAction
{
  public:
    ExN02EventAction();
    ~ExN02EventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();
    
  private:

};

#endif

    
