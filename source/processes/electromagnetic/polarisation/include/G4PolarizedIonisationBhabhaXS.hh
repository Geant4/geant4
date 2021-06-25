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
// Geant4 class file
//
// File name:     G4PolarizedIonisationBhabhaXS
//
// Author:        Andreas Schaelicke and Pavel Starovoitov
//
// Class Description:
//   * calculates the differential cross section (ME squared, without phase
//     space)
//   * incoming positron (along positive z direction) scatters at an electron at
//     rest
//   * phi denotes the angle between the scattering plane and the incoming PRF
//     X-axis
//   * all stokes vectors refer to spins in the Global System (X,Y,Z)

#ifndef G4PolarizedIonisationBhabhaXS_h
#define G4PolarizedIonisationBhabhaXS_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedXS.hh"

class G4PolarizedIonisationBhabhaXS : public G4VPolarizedXS
{
 public:
  G4PolarizedIonisationBhabhaXS();
  ~G4PolarizedIonisationBhabhaXS() override;

  void Initialize(G4double x, G4double y, G4double phi,
                  const G4StokesVector& p0, const G4StokesVector& p1,
                  G4int flag = 0) override;

  G4double XSection(const G4StokesVector& pol2,
                    const G4StokesVector& pol3) override;
  G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
                         const G4StokesVector& pol0,
                         const G4StokesVector& pol1) override;

  // return expected mean polarisation
  G4StokesVector GetPol2() override;
  G4StokesVector GetPol3() override;

  G4PolarizedIonisationBhabhaXS& operator                             =(
    const G4PolarizedIonisationBhabhaXS& right) = delete;
  G4PolarizedIonisationBhabhaXS(const G4PolarizedIonisationBhabhaXS&) = delete;

 private:
  // - part depending on the polarization of the final positron
  G4ThreeVector fPhi2;
  // - part depending on the polarization of the final electron
  G4ThreeVector fPhi3;

  G4double fPhi0;
};
#endif
