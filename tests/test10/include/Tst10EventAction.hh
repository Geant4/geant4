// $Id: Tst10EventAction.hh,v 1.2 1999-04-17 08:01:45 kurasige Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserEventAction
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Tst10EventAction_h
#define Tst10EventAction_h 1

#include "G4UserEventAction.hh"

class G4Event;

class Tst10EventAction : public G4UserEventAction
{
  public:
    Tst10EventAction();
    virtual ~Tst10EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
};

#endif

    
