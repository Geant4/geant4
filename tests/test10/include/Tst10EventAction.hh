// $Id: Tst10EventAction.hh,v 1.1 1999-01-08 16:35:32 gunter Exp $
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
    ~Tst10EventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();
};

#endif

    
