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
// $Id: G4ITTrackHolder.cc 64374 2012-10-31 16:37:23Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
//
// History:
// -----------
// 16 Mai 2012 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITTrackHolder.hh"
#include "G4Track.hh"

G4ITTrackHolder* G4ITTrackHolder::fInstance(0);

G4ITTrackHolder::G4ITTrackHolder()
{
    fInstance = this;
}

G4ITTrackHolder::~G4ITTrackHolder()
{
}

G4ITTrackHolder* G4ITTrackHolder::Instance()
{
    return fInstance;
}

void G4ITTrackHolder::PushTrack(G4Track* track)
{
    delete track;
}

double G4ITTrackHolder::GetTimeStep() const
{
    return 0;
}
