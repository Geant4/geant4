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
#ifndef Tst21EventAction_h
#define Tst21EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class Tst21EventAction : public G4UserEventAction
{
  public:
    Tst21EventAction(){evnum=0;}
    ~Tst21EventAction(){}

  public:
    virtual void BeginOfEventAction(const G4Event*)
    {
      G4cout <<"Event number "<<evnum++<<G4endl;
    }
    virtual void EndOfEventAction(const G4Event*){ }
 
  private:
    int evnum;
};

#endif

    
