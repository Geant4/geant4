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
// $Id: RE01SteppingAction.hh,v 1.1 2004/11/26 07:37:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//


#ifndef RE01SteppingAction_H
#define RE01SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class RE01SteppingAction : public G4UserSteppingAction
{
  public:
    RE01SteppingAction();
    virtual ~RE01SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

