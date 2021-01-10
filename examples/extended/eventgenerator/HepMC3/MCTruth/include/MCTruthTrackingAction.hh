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
/// \file eventgenerator/HepMC/MCTruth/include/MCTruthTrackingAction.hh
/// \brief Definition of the MCTruthTrackingAction class
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - MCTruthTrackingAction class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
// Date  : 2005-08-24
//
// --------------------------------------------------------------
#ifndef INCLUDE_MCTRUTHTRACKINGACTION_H 
#define INCLUDE_MCTRUTHTRACKINGACTION_H 1

#include "G4UserTrackingAction.hh"
#include "MCTruthManager.hh"

class MCTruthTrackingAction : public G4UserTrackingAction
{
public: 

  MCTruthTrackingAction( ); 

  virtual ~MCTruthTrackingAction( );
  
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

private:

  G4LorentzVector fmom;
  G4bool TrackToBeStored(const G4Track*);

};

#endif // INCLUDE_MCTRUTHTRACKINGACTION_H
