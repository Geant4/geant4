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
#ifndef HECEventAction_h
#define HECEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class HECEventAction : public G4UserEventAction
{
  public:
    HECEventAction();
    ~HECEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    void BadEvent(){GoodEv = false;}
 
  private:
    static int evnum;
    static G4bool GoodEv;
    static G4bool store;
    G4int moveCollID;
    G4int frontCollID;
    G4int larCollID;
    G4int cuCollID;
    G4int leakCollID;
};

#endif

    
