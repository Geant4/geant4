// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVASteppingAction.hh,v 1.3 2001-02-07 17:31:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVASteppingAction_h
#define TstVASteppingAction_h 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "Histo.hh"

class TstVASteppingAction : public G4UserSteppingAction
{
  public:
    TstVASteppingAction();
    virtual ~TstVASteppingAction();

    virtual void UserSteppingAction(const G4Step*);
private:
  odHisto Steplength;
  odHisto SteplengthProfile;
};

#endif
