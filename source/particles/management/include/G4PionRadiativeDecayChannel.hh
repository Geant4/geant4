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
// G4PionRadiativeDecayChannel
//
// Class decription:
//
// This class describes radiative pion decay kinematics.
// Rate is for gammas > 100keV.
// NOTE: it gives incorrect energy spectrum for neutrino

// Author: P.Gumplinger, 30 July 2007
// Reference: M. Blecher, TRIUMF/PIENU Technote
//            "Inclusion of pi->enug in the Monte Carlo"
// --------------------------------------------------------------------
#ifndef G4PionRadiativeDecayChannel_hh
#define G4PionRadiativeDecayChannel_hh 1

#include "G4ThreeVector.hh"
#include "G4VDecayChannel.hh"
#include "Randomize.hh"
#include "globals.hh"

class G4PionRadiativeDecayChannel : public G4VDecayChannel
{
  public:
    G4PionRadiativeDecayChannel(const G4String& theParentName, G4double theBR);
    ~G4PionRadiativeDecayChannel() override = default;

    G4DecayProducts* DecayIt(G4double) override;

  protected:
    G4PionRadiativeDecayChannel() = default;

    // Copy constructor and assignment operator
    G4PionRadiativeDecayChannel(const G4PionRadiativeDecayChannel&) = default;
    G4PionRadiativeDecayChannel& operator=(const G4PionRadiativeDecayChannel&);
};

#endif
