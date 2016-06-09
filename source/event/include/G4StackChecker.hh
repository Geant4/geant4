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
// $Id: G4StackChecker.hh,v 1.1 2003/06/12 15:14:43 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#ifndef G4StackChecker_h
#define G4StackChecker_h 1

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4ThreeVector.hh"

// class description:
//
//  This is the UserStackingAction for removing unreasonable tracks.
//

class G4Track;

class G4StackChecker  : public G4UserStackingAction
{
public:
  G4StackChecker();
  virtual ~G4StackChecker();

public: // with description

  virtual G4ClassificationOfNewTrack
        ClassifyNewTrack(const G4Track* track);

  virtual void NewStage() {};
  virtual void PrepareNewEvent() {};

private:

  G4ThreeVector nullDirection;

};

#endif
