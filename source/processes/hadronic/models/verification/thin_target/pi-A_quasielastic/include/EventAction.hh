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
// $Id: EventAction.hh,v 1.1 2003-07-31 01:15:59 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"


class RunAction;
class EventActionMessenger;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction*);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    
    void StoreDynamicParticle(const G4DynamicParticle*);
    void SetPrintModulo(G4int val)  {printModulo = val;};
    
  private:

    G4int printModulo;
    EventActionMessenger* eventMessenger;                         
    RunAction* theRunAction;
    G4std::vector<G4DynamicParticle*> particleList;
};

#endif


