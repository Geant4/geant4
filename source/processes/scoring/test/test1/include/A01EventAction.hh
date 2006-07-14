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
// $Id: A01EventAction.hh,v 1.1 2006-07-14 14:42:49 asaim Exp $
// --------------------------------------------------------------
//
#ifndef A01EventAction_h
#define A01EventAction_h 1


#include "G4UserEventAction.hh"
#include "G4THitsMap.hh"
#include "globals.hh"

class A01EventActionMessenger;

class A01EventAction : public G4UserEventAction
{
  public:
    A01EventAction(G4bool);
    virtual ~A01EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int realWorldID[4];
    G4int paraWorldID[4];
    G4bool ifPara;
  private:
    G4double Total(G4THitsMap<G4double>*);
};

#endif
