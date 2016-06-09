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
// $Id: G4ITTrackHolder.hh 64374 2012-10-31 16:37:23Z gcosmo $
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

#ifndef G4ITTRACKHOLDER_HH
#define G4ITTRACKHOLDER_HH

class G4Track;

/**
  * G4ITTrackHolder is an empty interface that permits to push tracks to the IT system
  * without actually depending on the IT tracking system.
  * However, G4ITTrackHolder does not permit to retrieve any track.
  */

class G4ITTrackHolder
{
protected :
    G4ITTrackHolder();
    virtual ~G4ITTrackHolder();
    static G4ITTrackHolder* fInstance;

public:
    static G4ITTrackHolder* Instance();
    virtual void PushTrack(G4Track*);
    virtual double GetTimeStep() const ;
};

#endif // G4ITTRACKHOLDER_HH
