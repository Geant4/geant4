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
//
/// \file Par02TrackingAction.hh
/// \brief Definition of the Par02TrackingAction class

#ifndef PAR02_TRACKING_ACTION_H
#define PAR02_TRACKING_ACTION_H

#include "G4UserTrackingAction.hh"
#include "globals.hh"

/// Tracking action (before/after track processing).
///
/// Defines the action at the start and at the end of processing of each track.
/// The class needs to be set in G4RunManager::SetUserAction().
/// @author Anna Zaborowska

class Par02TrackingAction : public G4UserTrackingAction {
  public:
    
    /// A default constructor.
    Par02TrackingAction();

    virtual ~Par02TrackingAction();

    /// Defines the actions at the start of processing the track.
    /// It checks the pseudorapidity range and if the particle is a primary.
    virtual void  PreUserTrackingAction( const G4Track* track );
    
    /// Defines the actions at the end of processing the track. 
    /// It saves the information of MC data (PDG code, initial momentum),
    /// tracker (momentum), EMCal and HCal (energy deposit and its position)
    /// as well as resolution and efficiency for all the detectors.
    virtual void  PostUserTrackingAction( const G4Track* track );
};

#endif

