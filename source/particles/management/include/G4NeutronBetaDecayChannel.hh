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
// G4NeutronBetaDecayChannel
//
// Class decription:
//
// This class describes free neutron beta decay kinematics.
// It neglects neutron/electron polarization without Coulomb effect.

// Author: H.Kurashige, 18 September 2001
// --------------------------------------------------------------------
#ifndef G4NeutronBetaDecayChannel_hh
#define G4NeutronBetaDecayChannel_hh 1

#include "G4VDecayChannel.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4NeutronBetaDecayChannel : public G4VDecayChannel
{
  public:
    G4NeutronBetaDecayChannel(const G4String& theParentName, G4double theBR);
    ~G4NeutronBetaDecayChannel() override = default;

    G4DecayProducts* DecayIt(G4double) override;

  protected:
    G4NeutronBetaDecayChannel() = default;

    // Copy constructor and assignment operator
    G4NeutronBetaDecayChannel(const G4NeutronBetaDecayChannel&);
    G4NeutronBetaDecayChannel& operator=(const G4NeutronBetaDecayChannel&);

    // e-neutrino angular correlation parameter
    const G4double aENuCorr = -0.102;
};

#endif
