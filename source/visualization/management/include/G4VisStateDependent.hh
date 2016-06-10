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
// $Id: G4VisStateDependent.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// A "state dependent" service class for G4VisManager.
// John Allison  29th November 1999

// G4VisStateDependent is "state dependent", i.e., it is notified on
// change of state (G4ApplicationState).  This is used to message the
// G4VisManager to draw hits and trajectories in the current scene at the
// end of event, as required.

#ifndef G4VISSTATEDEPENDENT_HH
#define G4VISSTATEDEPENDENT_HH

#include "G4VStateDependent.hh"

class G4VisManager;

class G4VisStateDependent: public G4VStateDependent {
  friend class G4VisManager;
private:
  G4VisStateDependent (G4VisManager *);
  G4bool Notify (G4ApplicationState requestedState);
  G4VisManager* fpVisManager;
};

#endif
