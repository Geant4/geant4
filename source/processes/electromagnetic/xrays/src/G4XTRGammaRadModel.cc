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

#include "G4XTRGammaRadModel.hh"

////////////////////////////////////////////////////////////////////////////
// Constructor, destructor
G4XTRGammaRadModel::G4XTRGammaRadModel(G4LogicalVolume* anEnvelope,
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
  fExitFlux = true;
}

///////////////////////////////////////////////////////////////////////////
G4XTRGammaRadModel::~G4XTRGammaRadModel() = default;

void G4XTRGammaRadModel::ProcessDescription(std::ostream& out) const
{
  out << "Rough model describing X-ray transition radiation. Thicknesses of "
         "plates\n"
         "and gas gaps are distributed according to gamma distributions.\n";
}

///////////////////////////////////////////////////////////////////////////
// Rough approximation for radiator interference factor for the case of
// fully GamDistr radiator. The plate and gas gap thicknesses are distributed
// according to exponent. The mean values of the plate and gas gap thicknesses
// are supposed to be about XTR formation zones but much less than
// mean absorption length of XTR photons in coresponding material.
G4double G4XTRGammaRadModel::GetStackFactor(G4double energy, G4double gamma,
                                            G4double varAngle)
{
  G4double result, Qa, Qb, Q, Za, Zb, Ma, Mb;

  Za = GetPlateFormationZone(energy, gamma, varAngle);
  Zb = GetGasFormationZone(energy, gamma, varAngle);

  Ma = GetPlateLinearPhotoAbs(energy);
  Mb = GetGasLinearPhotoAbs(energy);

  Qa = (1.0 + fPlateThick * Ma / fAlphaPlate);
  Qa = std::pow(Qa, -fAlphaPlate);
  Qb = (1.0 + fGasThick * Mb / fAlphaGas);
  Qb = std::pow(Qb, -fAlphaGas);
  Q  = Qa * Qb;

  G4complex Ca(1.0 + 0.5 * fPlateThick * Ma / fAlphaPlate,
               fPlateThick / Za / fAlphaPlate);
  G4complex Cb(1.0 + 0.5 * fGasThick * Mb / fAlphaGas,
               fGasThick / Zb / fAlphaGas);

  G4complex Ha = std::pow(Ca, -fAlphaPlate);
  G4complex Hb = std::pow(Cb, -fAlphaGas);
  G4complex H  = Ha * Hb;

  G4complex F1 = (0.5 * (1 + Qa) * (1.0 + H) - Ha - Qa * Hb) / (1.0 - H);

  G4complex F2 = (1.0 - Ha) * (Qa - Ha) * Hb / (1.0 - H) / (Q - H);

  F2 *= std::pow(Q, G4double(fPlateNumber)) - std::pow(H, fPlateNumber);

  result = (1. - std::pow(Q, G4double(fPlateNumber))) / (1. - Q);

  G4complex stack = result * F1;
  stack += F2;
  stack *= 2.0 * OneInterfaceXTRdEdx(energy, gamma, varAngle);

  result = std::real(stack);

  return result;
}
