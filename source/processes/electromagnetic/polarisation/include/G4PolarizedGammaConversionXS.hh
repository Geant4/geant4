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
// File name:     G4PolarizedGammaConversionXS
//
// Author:        Karim Laihem

#ifndef G4PolarizedGammaConversionXS_h
#define G4PolarizedGammaConversionXS_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedXS.hh"

class G4PolarizedGammaConversionXS : public G4VPolarizedXS
{
 public:
  G4PolarizedGammaConversionXS();
  ~G4PolarizedGammaConversionXS() override;

  void Initialize(G4double eps, G4double X, G4double phi,
                  const G4StokesVector& p0, const G4StokesVector& p1,
                  G4int flag = 0) override;

  G4double XSection(const G4StokesVector& pol2,
                    const G4StokesVector& pol3) override;

  // return expected mean polarisation
  G4StokesVector GetPol2() override;  // electron/positron
  G4StokesVector GetPol3() override;  // photon

  G4PolarizedGammaConversionXS& operator                            =(
    const G4PolarizedGammaConversionXS& right) = delete;
  G4PolarizedGammaConversionXS(const G4PolarizedGammaConversionXS&) = delete;

 private:
  static G4double SCRN[2][19];  // screening function lookup table

  G4StokesVector fFinalElectronPolarization;
  G4StokesVector fFinalPositronPolarization;
};

#endif
