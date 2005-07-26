//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

#ifndef G4HumanPhantomEventAction_h
#define G4HumanPhantomEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class G4HumanPhantomEventAction : public G4UserEventAction
{
  public:
    G4HumanPhantomEventAction();
   ~G4HumanPhantomEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    G4int GetEvent() const {return eventNumber;};
    void SetPath(G4double);

  private:
    G4int eventNumber;
    G4double path;
};
#endif

    
