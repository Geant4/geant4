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
// G4MuonDecayChannelWithSpin
//
// Class decription:
//
// This class describes muon decay kinematics.
// It assumes V-A coupling with 1st order radiative corrections, the
// standard model parameter values, but gives incorrect energy spectrum
// for neutrinos.
// References:
// - Florian Scheck "Muon Physics", in Physics Reports
//   (Review Section of Physics Letters) 44, No. 4 (1978)
//   187-248. North-Holland Publishing Company, Amsterdam at page 210 cc.
// - W.E. Fisher and F. Scheck, Nucl. Phys. B83 (1974) 25.

// Authors: P.Gumplinger and T.MacPhail, 17 August 2004
// --------------------------------------------------------------------
#ifndef G4MuonDecayChannelWithSpin_hh
#define G4MuonDecayChannelWithSpin_hh 1

#include "G4MuonDecayChannel.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <CLHEP/Units/PhysicalConstants.h>

class G4MuonDecayChannelWithSpin : public G4MuonDecayChannel
{
  public:
    G4MuonDecayChannelWithSpin(const G4String& theParentName, G4double theBR);
    ~G4MuonDecayChannelWithSpin() override = default;

    G4DecayProducts* DecayIt(G4double) override;

  protected:
    // Copy constructor and assignment operator
    G4MuonDecayChannelWithSpin(const G4MuonDecayChannelWithSpin&) = default;
    G4MuonDecayChannelWithSpin& operator=(const G4MuonDecayChannelWithSpin&);

  private:
    G4MuonDecayChannelWithSpin() = default;

    // Radiative Correction Factors
    inline G4double F_c(G4double x, G4double x0, G4double omega);
    inline G4double F_theta(G4double x, G4double x0, G4double omega);
    G4double R_c(G4double x, G4double omega);
};

// ------------------------
// Inline methods
// ------------------------

inline G4double G4MuonDecayChannelWithSpin::F_c(G4double x, G4double x0, G4double omega)
{
  G4double f_c;

  f_c = (5. + 17. * x - 34. * x * x) * (omega + std::log(x)) - 22. * x + 34. * x * x;
  f_c = (1. - x) / (3. * x * x) * f_c;
  f_c = (6. - 4. * x) * R_c(x, omega) + (6. - 6. * x) * std::log(x) + f_c;
  f_c = (CLHEP::fine_structure_const / CLHEP::twopi) * (x * x - x0 * x0) * f_c;

  return f_c;
}

inline G4double G4MuonDecayChannelWithSpin::F_theta(G4double x, G4double x0, G4double omega)
{
  G4double f_theta;

  f_theta = (1. + x + 34 * x * x) * (omega + std::log(x)) + 3. - 7. * x - 32. * x * x;
  f_theta = f_theta + ((4. * (1. - x) * (1. - x)) / x) * std::log(1. - x);
  f_theta = (1. - x) / (3. * x * x) * f_theta;
  f_theta = (2. - 4. * x) * R_c(x, omega) + (2. - 6. * x) * std::log(x) - f_theta;
  f_theta = (CLHEP::fine_structure_const / CLHEP::twopi) * (x * x - x0 * x0) * f_theta;

  return f_theta;
}

#endif
