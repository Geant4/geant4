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
/*
 * G4ITLeadingTracks.cc
 *
 *  Created on: Jun 26, 2015
 *      Author: mkaramit
 */

#include <G4ITLeadingTracks.hh>
#include "G4TrackingInformation.hh"
#include <G4IT.hh>

G4ITLeadingTracks::G4ITLeadingTracks()
{
  // TODO Auto-generated constructor stub

}

G4ITLeadingTracks::~G4ITLeadingTracks()
{
  // TODO Auto-generated destructor stub
}


void G4ITLeadingTracks::Reset()
{
  if (fLeadingTracks.empty() == false)
  {
    std::vector<G4Track*>::iterator fLeadingTracks_i = fLeadingTracks.begin();

    while (fLeadingTracks_i != fLeadingTracks.end())
    {
      G4Track* track = *fLeadingTracks_i;
      if (track)
      {
        G4IT* ITrack = GetIT(track);
        if (ITrack)
        {
          ITrack->GetTrackingInfo()->SetLeadingStep(false);
        }
      }

      ++fLeadingTracks_i;
      continue;
    }

    fLeadingTracks.clear();
  }
}

void G4ITLeadingTracks::Push(G4Track* track)
{
  fLeadingTracks.push_back(track);
}

void G4ITLeadingTracks::PrepareLeadingTracks()
{
  for(size_t i = 0 ; i < fLeadingTracks.size() ; ++i)
  {
    G4Track* track = fLeadingTracks[i];
    G4IT* ITrack = GetIT(track);
    ITrack->GetTrackingInfo()->SetLeadingStep(true);
//    ITrack->GetTrackingInfo()->SetLeadingStep(false);
  }
}
