// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingAction.hh,v 2.1 1998/07/12 02:37:19 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
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
