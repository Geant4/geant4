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
// $Id: PersEx02EventAction.hh,v 1.4 2001/07/11 09:58:14 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
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

    
