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

// ====================================================================
//
//   XXEventAction.hh
//   $Id: XXEventAction.hh,v 1.1 2002-04-29 20:44:34 asaim Exp $
//
// ====================================================================
#ifndef XX_EVENT_ACTION_H
#define XX_EVENT_ACTION_H

#include "G4UserEventAction.hh"

class XXEventAction : public G4UserEventAction {
public:
  XXEventAction();
  ~XXEventAction();

  virtual void BeginOfEventAction(const G4Event* anEvent);
  virtual void EndOfEventAction(const G4Event* anEvent);
};

#endif
