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
// $ID:      Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef Tst32SteppingAction_H
#define Tst32SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst32SteppingAction : public G4UserSteppingAction
{
  public:
    Tst32SteppingAction();
    ~Tst32SteppingAction();
    virtual void UserSteppingAction(const G4Step*);
private:
  G4int nevent;
};

#endif

