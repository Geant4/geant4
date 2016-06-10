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
// $Id: G4ITTrackingInteractivity.cc 91584 2015-07-27 13:01:48Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
//
// History:
// -----------
//
//
// -------------------------------------------------------------------

#include <G4ITSteppingVerbose.hh>
#include "G4ITTrackingInteractivity.hh"
#include "G4Track.hh"
#include "G4String.hh"

void G4ITTrackingInteractivity::TrackBanner(G4Track* track, const G4String& message)
{
  G4cout << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  if(message != "")
      G4cout << message ;
  G4cout << " * G4Track Information: "
         << "   Particle : " << track->GetDefinition()->GetParticleName()
         << ","
         << "   Track ID : " << track->GetTrackID()
         << ","
         << "   Parent ID : " << track->GetParentID()
         << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << G4endl;
}

G4ITTrackingInteractivity::G4ITTrackingInteractivity(G4VITSteppingVerbose* verbose)
{
  fVerboseLevel = 0;
  if(verbose)
  {
    fpVerbose = verbose;
  }
  else
  {
    fpVerbose = new G4ITSteppingVerbose();
  }
}

void G4ITTrackingInteractivity::SetSteppingVerboseLevel(G4int level)
{
  fpVerbose->SetVerbose(level);
}

//void G4ITTrackingInteractivity::StartTracking(G4Track* track)
//{
//  fpVerbose->TrackingStarted(track);
//}
//
//void G4ITTrackingInteractivity::EndTracking(G4Track* track)
//{
//  fpVerbose->TrackingEnded(track);
//}

G4int G4ITTrackingInteractivity::GetSteppingVerboseLevel() const
{
  return fpVerbose->GetVerbose();
}
