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

#ifndef G4VUPPACTION_H
#define G4VUPPACTION_H


#include "globals.hh"
#include "G4UppInteraction.hh"
#include "G4UppTrackVector.hh"


class G4VUppAction
{
public:

  virtual G4UppTrackChange* perform(const G4UppTrackVector& tracks) const = 0;

  virtual G4bool isValid() const = 0;

  virtual void setActionTime(const G4double newTime) 
    { actionTime = newTime; }
  virtual G4double getActionTime() const 
    { return actionTime; }

  virtual void dump() const 
    { G4cout << "Unknown Action" << G4endl; }

private:

  G4double actionTime;

};


#endif // G4VUPPACTION_H
