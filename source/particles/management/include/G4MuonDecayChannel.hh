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
// G4MuonDecayChannel
//
// Class decription:
//
// This class describes muon decay kinematics.
// This version neglects muon polarization; assumes the pure
// V-A coupling gives incorrect energy spectrum for neutrinos.

// Author: H.Kurashige, 30 May 1997
// --------------------------------------------------------------------
#ifndef G4MuonDecayChannel_hh
#define G4MuonDecayChannel_hh 1

#include "G4VDecayChannel.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4MuonDecayChannel : public G4VDecayChannel
{
  public:
    // Constructor & destructor
    G4MuonDecayChannel(const G4String& parentName, G4double BR);
    ~G4MuonDecayChannel() override = default;

    G4DecayProducts* DecayIt(G4double) override;

  protected:
    G4MuonDecayChannel() = default;

    // Copy constructor and assignment operator
    G4MuonDecayChannel(const G4MuonDecayChannel&) = default;
    G4MuonDecayChannel& operator=(const G4MuonDecayChannel&);
};

#endif
