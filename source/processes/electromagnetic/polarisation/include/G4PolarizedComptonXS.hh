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
// Geant4 Class file
//
// File name:     G4PolarizedComptonXS
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   determine the  polarization of the final state
//   in a Compton scattering process employing the differential
//   cross section by F.W.Lipps & H.A.Tolhoek
//   ( Physica 20 (1954) 395 )

#ifndef G4PolarizedComptonXS_h
#define G4PolarizedComptonXS_h 1

#include "globals.hh"
#include "G4StokesVector.hh"
#include "G4VPolarizedXS.hh"

class G4PolarizedComptonXS : public G4VPolarizedXS
{
 public:
  G4PolarizedComptonXS();
  ~G4PolarizedComptonXS() override;

  // prepares the ingredients for the calculation of a polarization
  // dependent differential cross section
  // the kinematics is fixed (X - incoming photon energy in units of electron
  // mass, eps - outgoing photon energy in unit of incoming photon energy,
  // and polarization of the incoming particles fixed (p0, p1)
  // a flag specifies the extent to which polarization is taken into account
  void Initialize(G4double eps, G4double X, G4double phi,
                  const G4StokesVector& p0, const G4StokesVector& p1,
                  G4int flag = 0) override;

  // returns the differential cross section for a given polarisation state
  // of the final state particles to be used in the calculation of the
  // polarization transfer the calculation has to be initialised by calling
  // Initialize() prior to the first call of this function (see above)
  G4double XSection(const G4StokesVector& pol2,
                    const G4StokesVector& pol3) override;
  // total cross section
  G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
                         const G4StokesVector& pol0,
                         const G4StokesVector& pol1) override;

  // return expected mean polarisation
  G4StokesVector GetPol2() override;
  G4StokesVector GetPol3() override;

  G4PolarizedComptonXS& operator=(const G4PolarizedComptonXS& right) = delete;
  G4PolarizedComptonXS(const G4PolarizedComptonXS&)                  = delete;

 private:
  void DefineCoefficients(const G4StokesVector& pol0,
                          const G4StokesVector& pol1);

  static constexpr G4double re2 =
    CLHEP::classic_electr_radius * CLHEP::classic_electr_radius *
    (4. * CLHEP::pi / CLHEP::hbarc) * (4. * CLHEP::pi / CLHEP::hbarc);

  // these variables store the information necessary to evaluate the
  // differential cross section for arbitrary final state
  // polarizations (used in XSection):
  // - part depending on the polarization of the final photon
  G4ThreeVector fPhi2;
  // - part depending on the polarization of the final electron
  G4ThreeVector fPhi3;
  // - polarization independent part
  G4double fPhi0;
  // - product of polarizations of initial particles
  G4double polxx, polyy, polzz, polxz, polzx, polyz, polzy, polxy, polyx;

  // G4double diffXSFactor, totalXSFactor;
  G4double fPolXS, fUnpXS;
};

#endif
