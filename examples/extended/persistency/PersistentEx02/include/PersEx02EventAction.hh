// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02EventAction.hh,v 1.3 1999/11/29 18:23:31 morita Exp $
// GEANT4 tag $Name: geant4-01-00 $
//


#ifndef PersEx02EventAction_h
#define PersEx02EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class PersEx02EventAction : public G4UserEventAction
{
  public:
    PersEx02EventAction();
    virtual ~PersEx02EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int colID1;
    G4int colID2;
};

#endif

    
