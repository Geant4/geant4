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

#include "G4XTRRegularRadModel.hh"

#include "G4PhysicalConstants.hh"

////////////////////////////////////////////////////////////////////////////
// Constructor, destructor
G4XTRRegularRadModel::G4XTRRegularRadModel(G4LogicalVolume* anEnvelope,
                                           G4Material* foilMat,
                                           G4Material* gasMat, G4double a,
                                           G4double b, G4int n,
                                           const G4String& processName)
  : G4VXTRenergyLoss(anEnvelope, foilMat, gasMat, a, b, n, processName)
{
  G4cout << " XTR Regular discrete radiator model is called" << G4endl;

  fExitFlux = true;
}

///////////////////////////////////////////////////////////////////////////
G4XTRRegularRadModel::~G4XTRRegularRadModel() = default;

///////////////////////////////////////////////////////////////////////////
void G4XTRRegularRadModel::ProcessDescription(std::ostream& out) const
{
  out << "Describes X-ray transition radiation with thickness of gaps and "
         "plates\n"
         "fixed.\n";
}

///////////////////////////////////////////////////////////////////////////
G4double G4XTRRegularRadModel::SpectralXTRdEdx(G4double energy)
{
  static constexpr G4double cofPHC = 4. * pi * hbarc;
  G4double result, sum = 0., tmp, cof1, cof2, cofMin, theta2, theta2k;
  G4double aMa, bMb, sigma, dump;
  G4int k, kMax, kMin;

  aMa   = fPlateThick * GetPlateLinearPhotoAbs(energy);
  bMb   = fGasThick * GetGasLinearPhotoAbs(energy);
  sigma = 0.5 * (aMa + bMb);
  dump  = std::exp(-fPlateNumber * sigma);
  if(verboseLevel > 2)
    G4cout << " dump = " << dump << G4endl;
  tmp  = (fSigma1 - fSigma2) / cofPHC / energy;
  cof1 = fPlateThick * tmp;
  cof2 = fGasThick * tmp;

  cofMin = energy * (fPlateThick + fGasThick) / fGamma / fGamma;
  cofMin += (fPlateThick * fSigma1 + fGasThick * fSigma2) / energy;
  cofMin /= cofPHC;

  theta2 = cofPHC / (energy * (fPlateThick + fGasThick));

  kMin = G4int(cofMin);
  if(cofMin > kMin)
    kMin++;

  kMax = kMin + 49;

  if(verboseLevel > 2)
  {
    G4cout << cof1 << "     " << cof2 << "        " << cofMin << G4endl;
    G4cout << "kMin = " << kMin << ";    kMax = " << kMax << G4endl;
  }
  for(k = kMin; k <= kMax; ++k)
  {
    tmp    = pi * fPlateThick * (k + cof2) / (fPlateThick + fGasThick);
    result = (k - cof1) * (k - cof1) * (k + cof2) * (k + cof2);
    if(k == kMin && kMin == G4int(cofMin))
    {
      sum +=
        0.5 * std::sin(tmp) * std::sin(tmp) * std::abs(k - cofMin) / result;
    }
    else
    {
      sum += std::sin(tmp) * std::sin(tmp) * std::abs(k - cofMin) / result;
    }
    theta2k = std::sqrt(theta2 * std::abs(k - cofMin));

    if(verboseLevel > 2)
    {
      G4cout << k << "   " << theta2k << "     "
             << std::sin(tmp) * std::sin(tmp) * std::abs(k - cofMin) / result
             << "      " << sum << G4endl;
    }
  }
  result = 2 * (cof1 + cof2) * (cof1 + cof2) * sum / energy;
  result *= dump * (-1 + dump + 2 * fPlateNumber);

  return result;
}

///////////////////////////////////////////////////////////////////////////
// Approximation for radiator interference factor for the case of
// fully Regular radiator. The plate and gas gap thicknesses are fixed.
// The mean values of the plate and gas gap thicknesses
// are supposed to be about XTR formation zones but much less than
// mean absorption length of XTR photons in corresponding material.
G4double G4XTRRegularRadModel::GetStackFactor(G4double energy, G4double gamma,
                                              G4double varAngle)
{
  G4double aZa = fPlateThick / GetPlateFormationZone(energy, gamma, varAngle);
  G4double bZb = fGasThick / GetGasFormationZone(energy, gamma, varAngle);

  G4double aMa = fPlateThick * GetPlateLinearPhotoAbs(energy);
  G4double bMb = fGasThick * GetGasLinearPhotoAbs(energy);

  G4double Qa = std::exp(-aMa);
  G4double Qb = std::exp(-bMb);
  G4double Q  = Qa * Qb;

  G4complex Ha(std::exp(-0.5 * aMa) * std::cos(aZa),
               -std::exp(-0.5 * aMa) * std::sin(aZa));

  G4complex Hb(std::exp(-0.5 * bMb) * std::cos(bZb),
               -std::exp(-0.5 * bMb) * std::sin(bZb));

  G4complex H  = Ha * Hb;
  G4complex Hs = std::conj(H);

  G4complex F2 = (1.0 - Ha) * (Qa - Ha) * Hb * (1.0 - Hs) * (Q - Hs);
  F2 *= std::pow(Q, G4double(fPlateNumber)) - std::pow(H, fPlateNumber);

  G4double result = (1. - std::pow(Q, G4double(fPlateNumber))) / (1. - Q);
  result *= (1. - Qa) * (1. + Qa - 2. * std::sqrt(Qa) * std::cos(aZa));
  result /= (1. - std::sqrt(Q)) * (1. - std::sqrt(Q)) +
            4. * std::sqrt(Q) * std::sin(0.5 * (aZa + bZb)) *
              std::sin(0.5 * (aZa + bZb));

  G4double I2 = 1.;
  I2 /= (1. - std::sqrt(Q)) * (1. - std::sqrt(Q)) +
        4. * std::sqrt(Q) * std::sin(0.5 * (aZa + bZb)) *
          std::sin(0.5 * (aZa + bZb));

  I2 /= Q * ((std::sqrt(Q) - std::cos(aZa + bZb)) *
               (std::sqrt(Q) - std::cos(aZa + bZb)) +
             std::sin(aZa + bZb) * std::sin(aZa + bZb));

  G4complex stack = 2. * I2 * F2;
  stack += result;
  stack *= OneInterfaceXTRdEdx(energy, gamma, varAngle);

  return std::real(stack);
}
