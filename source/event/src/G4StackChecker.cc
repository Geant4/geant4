//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4StackChecker
//
// Author: Makoto Asai, 2003
// --------------------------------------------------------------------
#include "G4StackChecker.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4ios.hh"

G4ClassificationOfNewTrack
G4StackChecker::ClassifyNewTrack(const G4Track* track)
{
  G4ClassificationOfNewTrack result = fUrgent;
  G4double e = track->GetKineticEnergy();
  if ( (!(e < 0.0) && !(e > 0.0) && !(e == 0.0))
    || track->GetMomentumDirection() == nullDirection)
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
