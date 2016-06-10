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
// $Id: G4VStateDependent.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//
//      ---------------- G4VStateDependent ----------------
//
// Authors: G.Cosmo, M.Asai - November 1996
//
// ------------------------------------------------------------
//
// Class description:
//
// Abstract base class of all classes which need to be notified when
// the state of Geant4 changes. The concrete class object derived from
// this base class will be automatically registered to G4StateManager
// and the virtual method Notify() will be invoked when the state changes.

// ------------------------------------------------------------

#ifndef G4VStateDependent_h
#define G4VStateDependent_h 1

#include "globals.hh"
#include "G4ApplicationState.hh"

class G4VStateDependent
{

public:

  explicit G4VStateDependent(G4bool bottom=false);
  virtual ~G4VStateDependent();
  G4int operator==(const G4VStateDependent &right) const;
  G4int operator!=(const G4VStateDependent &right) const;

public: // with description

  virtual G4bool Notify(G4ApplicationState requestedState) = 0;
    // Pure virtual method which will be invoked by G4StateManager.
    // In case state change must not be allowed by some reason of the
    // concrete class, false should be returned. But this scheme is
    // NOT recommended to use. All command which are state sensitive
    // MUST assign available state(s).

private:

  G4VStateDependent(const G4VStateDependent &right);
  G4VStateDependent& operator=(const G4VStateDependent &right);

};

#endif
