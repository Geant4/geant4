// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyEventAction.hh,v 1.1 1999-01-07 16:05:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyEventAction_h
#define MyEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction();
    ~MyEventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();

  private:
    G4int calorimeterCollID;
    /////////G4int bulkCollID;
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif

    
