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
// $Id: G4StackChecker.cc,v 1.1 2003/06/12 15:14:44 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//

#include "G4StackChecker.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4ios.hh"

G4StackChecker::G4StackChecker():
   nullDirection(G4ThreeVector(0.0,0.0,0.0))
{}

G4StackChecker::~G4StackChecker()
{}

G4ClassificationOfNewTrack G4StackChecker::ClassifyNewTrack
(const G4Track* track)
{
  G4ClassificationOfNewTrack result = fUrgent;
  G4double e = track->GetKineticEnergy();
  if ( (!(e < 0.0) && !(e > 0.0) && !(e == 0.0)) ||
       track->GetMomentumDirection() == nullDirection)
    {
       result = fKill;
       G4String nam = track->GetDefinition()->GetParticleName();
       G4cout << "### G4StackChecker: event# "
              << (G4EventManager::GetEventManager())->GetConstCurrentEvent()->GetEventID()
	      << " unacceptable " << nam << " is killed in the stack" << G4endl;
       G4cout << "### " << nam << " have been produced by the process " 
              <<  track->GetCreatorProcess()->GetProcessName() 
              << " trackID= " << track->GetTrackID()
              << " parentID= " << track->GetParentID()
              << G4endl;
       G4cout << "### E= " << track->GetKineticEnergy()
	      << " position= " << track->GetPosition()
	      << " direction= " << track->GetMomentumDirection()
	      << " time= " << track->GetGlobalTime()
	      << G4endl;
    }
  return result;
}

