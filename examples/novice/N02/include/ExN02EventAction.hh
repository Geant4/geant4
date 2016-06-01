// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02EventAction.hh,v 1.3 1998/10/09 14:31:02 japost Exp $
// GEANT4 tag $Name: geant4-00 $
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

    
