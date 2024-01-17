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
// G4TauLeptonicDecayChannel
//
// Class decription:
//
// Class describing tau leptonic decay kinematics.
// This version assumes the pure V-A coupling and gives incorrect
// energy spectrum for neutrinos without tau polarization.

// Author: H.Kurashige, 30 May 1997
// --------------------------------------------------------------------
#ifndef G4TauLeptonicDecayChannel_hh
#define G4TauLeptonicDecayChannel_hh 1

#include "G4VDecayChannel.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4TauLeptonicDecayChannel : public G4VDecayChannel
{
  public:
    G4TauLeptonicDecayChannel(const G4String& theParentName, G4double theBR,
                              const G4String& theLeptonName);
    ~G4TauLeptonicDecayChannel() override = default;

    G4DecayProducts* DecayIt(G4double) override;

  protected:
    G4TauLeptonicDecayChannel() = default;

    G4TauLeptonicDecayChannel(const G4TauLeptonicDecayChannel&) = default;
    G4TauLeptonicDecayChannel& operator=(const G4TauLeptonicDecayChannel&);

  private:
    static G4double spectrum(G4double momentum, G4double energy, G4double mtau, G4double ml);
};

#endif
