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
// File name:     G4PolarizedPhotoElectricXS
//
// Author:        Andreas Schaelicke
//

#ifndef G4PolarizedPhotoElectricXS_h
#define G4PolarizedPhotoElectricXS_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedXS.hh"

class G4PolarizedPhotoElectricXS : public G4VPolarizedXS
{
 public:
  G4PolarizedPhotoElectricXS();
  ~G4PolarizedPhotoElectricXS() override;

  void Initialize(G4double aGammaE, G4double aLept0E, G4double sintheta,
                  const G4StokesVector& beamPol, const G4StokesVector&,
                  G4int flag = 0) override;

  G4double XSection(const G4StokesVector& pol2,
                    const G4StokesVector& pol3) override;

  // return expected mean polarisation
  G4StokesVector GetPol2() override;
  G4StokesVector GetPol3() override;

  G4PolarizedPhotoElectricXS& operator=(
    const G4PolarizedPhotoElectricXS& right) = delete;
  G4PolarizedPhotoElectricXS(const G4PolarizedPhotoElectricXS&) = delete;

 private:
  G4StokesVector fFinalElectronPolarization;
};

#endif
