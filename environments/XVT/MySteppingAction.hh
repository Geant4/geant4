// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingAction.hh,v 1.1 1999-01-07 16:05:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MySteppingAction_h
#define MySteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class MySteppingAction : public G4UserSteppingAction
{
  public:
    MySteppingAction();
    ~MySteppingAction(){};

    void UserSteppingAction();

  private:
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };

};

#endif
