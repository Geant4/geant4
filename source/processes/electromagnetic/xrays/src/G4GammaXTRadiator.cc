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

#include "G4GammaXTRadiator.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
// Constructor, destructor
G4GammaXTRadiator::G4GammaXTRadiator(G4LogicalVolume* anEnvelope,
                                     G4double alphaPlate, G4double alphaGas,
                                     G4Material* foilMat, G4Material* gasMat,
                                     G4double a, G4double b, G4int n,
                                     const G4String& processName)
  : G4VXTRenergyLoss(anEnvelope, foilMat, gasMat, a, b, n, processName)
{
  G4cout << "Gamma distributed X-ray TR radiator model is called" << G4endl;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  fAlphaPlate = alphaPlate;
  fAlphaGas   = alphaGas;
  G4cout << "fAlphaPlate = " << fAlphaPlate << " ; fAlphaGas = " << fAlphaGas
         << G4endl;
}

///////////////////////////////////////////////////////////////////////////
G4GammaXTRadiator::~G4GammaXTRadiator() = default;

void G4GammaXTRadiator::ProcessDescription(std::ostream& out) const
{
  out
    << "Rough approximation describing a radiator of X-ray transition "
       "radiation.\n"
       "Thicknesses of plates and gas gaps are distributed according to gamma\n"
       "description.\n";
}

///////////////////////////////////////////////////////////////////////////
// Rough approximation for radiator interference factor for the case of
// fully GamDistr radiator. The plate and gas gap thicknesses are distributed
// according to exponent. The mean values of the plate and gas gap thicknesses
// are supposed to be about XTR formation zones but much less than
// mean absorption length of XTR photons in corresponding material.
G4double G4GammaXTRadiator::GetStackFactor(G4double energy, G4double gamma,
                                           G4double varAngle)
{
  G4double result, Za, Zb, Ma, Mb;

  Za = GetPlateFormationZone(energy, gamma, varAngle);
  Zb = GetGasFormationZone(energy, gamma, varAngle);

  Ma = GetPlateLinearPhotoAbs(energy);
  Mb = GetGasLinearPhotoAbs(energy);

  G4complex Ca(1.0 + 0.5 * fPlateThick * Ma / fAlphaPlate,
               fPlateThick / Za / fAlphaPlate);
  G4complex Cb(1.0 + 0.5 * fGasThick * Mb / fAlphaGas,
               fGasThick / Zb / fAlphaGas);

  G4complex Ha = std::pow(Ca, -fAlphaPlate);
  G4complex Hb = std::pow(Cb, -fAlphaGas);
  G4complex H  = Ha * Hb;

  G4complex F1 = (1.0 - Ha) * (1.0 - Hb) / (1.0 - H) * G4double(fPlateNumber);

  G4complex F2 = (1.0 - Ha) * (1.0 - Ha) * Hb / (1.0 - H) / (1.0 - H) *
                 (1.0 - std::pow(H, fPlateNumber));

  G4complex R = (F1 + F2) * OneInterfaceXTRdEdx(energy, gamma, varAngle);

  result = 2.0 * std::real(R);

  return result;
}
