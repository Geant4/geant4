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
#ifndef Tst25StackingAction_h
#define Tst25StackingAction_h

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"

class Tst25StackingAction : public G4UserStackingAction
{
      virtual G4ClassificationOfNewTrack
        ClassifyNewTrack(const G4Track* aTrack)
	{
	  G4ClassificationOfNewTrack result(fUrgent);
// 	  if(aTrack->GetDefinition()->GetPDGCharge() == 0 && 
// 	    aTrack->GetDefinition()->GetPDGMass()<200*MeV)
// 	  {
// 	    result = fKill;
// 	  }
// 	  if(aTrack->GetKineticEnergy()<1*GeV)
// 	  {
// 	    result = fKill;
// 	  }
 	  return result;
	}
};

#endif
