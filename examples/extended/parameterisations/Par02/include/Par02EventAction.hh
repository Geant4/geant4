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
/// \file Par02EventAction.hh
/// \brief Definition of the Par02EventAction class

#ifndef PAR02_EVENT_ACTION_H
#define PAR02_EVENT_ACTION_H

#include "G4UserEventAction.hh"
#include "globals.hh"

/// Event action (before/after event processing).
///
/// Defines the action at the beginning and at the end of each event.
/// It is invoked by a G4EventManager when a G4Event object is sent
/// (which contains primary vertices and particles created by the 
/// Par02PrimaryGeneratorAction).
/// @author Anna Zaborowska

class Par02EventAction : public G4UserEventAction {
  public:
    
    /// A default constructor. 
    /// Sets the flag fSmear to true indicating that smearing will be performed.
    Par02EventAction();
    
    /// A constructor.
    /// @param aSmear The flag indicating if smearing has to be done.
    Par02EventAction( G4bool aSmear );

    virtual ~Par02EventAction();

    /// Defines the actions at the beginning of the event. 
    /// It sets the Par02EventInformation with fSmear flag. 
    /// It creates all the ntuples defined in Par02Output singleton class.
    virtual void BeginOfEventAction( const G4Event* );
    
    /// Defines the actions at the end of the event.
    virtual void EndOfEventAction( const G4Event* );

  private:
    
    /// A flag indicating if smearing should be performed. 
    /// Passed to Par02EventInformation in BeginOfEventAction(const G4Event*).
    G4bool fSmear;
};

#endif

