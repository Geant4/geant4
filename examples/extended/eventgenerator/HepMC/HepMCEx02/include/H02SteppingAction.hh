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
//   H02SteppingAction.hh
//   $Id: H02SteppingAction.hh,v 1.1 2002-05-28 14:10:53 murakami Exp $
//
// ====================================================================
#ifndef H02_STEPPING_ACTION_H
#define H02_STEPPING_ACTION_H

#include "G4UserSteppingAction.hh"

class H02SteppingAction : public G4UserSteppingAction {
public:
  H02SteppingAction();
  virtual ~H02SteppingAction();
  
  virtual void UserSteppingAction(const G4Step* astep);
};

#endif
