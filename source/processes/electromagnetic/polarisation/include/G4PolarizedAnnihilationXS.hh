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
// -------------------------------------------------------------------
//
// Geant4 Class file
//
// File name:     G4PolarizedAnnihilationXS
//
// Author:        Andreas Schaelicke and Pavel Starovoitov
//
// Class Description:
//   * calculates the differential cross section (ME squared,
//     without phase space) incoming positron (along positive z direction)
//     annihilations with an electron at rest
//   * phi denotes the angle between the scattering plane and
//     X axis of incoming partice reference frame (PRF)

#ifndef G4PolarizedAnnihilationXS_h
#define G4PolarizedAnnihilationXS_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedXS.hh"

class G4PolarizedAnnihilationXS : public G4VPolarizedXS
{
 public:
  G4PolarizedAnnihilationXS();
  virtual ~G4PolarizedAnnihilationXS() override;

  virtual void Initialize(G4double eps, G4double gamma, G4double phi,
                          const G4StokesVector& p0, const G4StokesVector& p1,
                          G4int flag = 0) override;

  G4double DiceEpsilon();

  virtual G4double XSection(const G4StokesVector& pol2,
                            const G4StokesVector& pol3) override;

  virtual G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
                                 const G4StokesVector& pol0,
                                 const G4StokesVector& pol1) override;

  // return expected mean polarisation
  virtual G4StokesVector GetPol2() override;
  virtual G4StokesVector GetPol3() override;

  // minimal energy fraction in TotalXSection
  virtual G4double GetXmin(G4double y) override;
  // maximal energy fraction in TotalXSection
  virtual G4double GetXmax(G4double y) override;

  G4double getVar(G4int);
  // test routine
  void getCoeff();

  G4PolarizedAnnihilationXS& operator=(const G4PolarizedAnnihilationXS& right) =
    delete;
  G4PolarizedAnnihilationXS(const G4PolarizedAnnihilationXS&) = delete;

 private:
  void TotalXS();
  void DefineCoefficients(const G4StokesVector& pol0,
                          const G4StokesVector& pol1);

  static constexpr G4double re2 =
    CLHEP::classic_electr_radius * CLHEP::classic_electr_radius;

  // - part depending on the polarization of the final positron
  G4ThreeVector fPhi2;

  // - part depending on the polarization of the final electron
  G4ThreeVector fPhi3;

  G4double polxx, polyy, polzz, polxz, polzx, polxy, polyx, polyz, polzy;

  // - unpolarised + part depending on the polarization of the initial pair
  G4double fPhi0;

  G4double fDice;
  G4double fPolXS, fUnpXS;
  G4double ISPxx, ISPyy, ISPzz, ISPnd;
};
#endif
