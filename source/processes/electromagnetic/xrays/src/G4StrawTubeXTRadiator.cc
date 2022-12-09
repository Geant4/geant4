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

#include "G4StrawTubeXTRadiator.hh"

#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

////////////////////////////////////////////////////////////////////////////
// Constructor, destructor
G4StrawTubeXTRadiator::G4StrawTubeXTRadiator(G4LogicalVolume* anEnvelope,
                                             G4Material* foilMat,
                                             G4Material* gasMat, G4double a,
                                             G4double b, G4Material* mediumMat,
                                             G4bool unishut,
                                             const G4String& processName)
  : G4VXTRenergyLoss(anEnvelope, foilMat, gasMat, a, b, 1, processName)
{
  if(verboseLevel > 0)
    G4cout << "Straw tube X-ray TR  radiator EM process is called" << G4endl;

  if(unishut)
  {
    fAlphaPlate = 1. / 3.;
    fAlphaGas   = 12.4;
    if(verboseLevel > 0)
      G4cout << "straw uniform shooting: "
             << "fAlphaPlate = " << fAlphaPlate
             << " ; fAlphaGas = " << fAlphaGas << G4endl;
  }
  else
  {
    fAlphaPlate = 0.5;
    fAlphaGas   = 5.;
    if(verboseLevel > 0)
      G4cout << "straw isotropical shooting: "
             << "fAlphaPlate = " << fAlphaPlate
             << " ; fAlphaGas = " << fAlphaGas << G4endl;
  }

  // index of medium material
  fMatIndex3 = (G4int)mediumMat->GetIndex();
  if(verboseLevel > 0)
    G4cout << "medium material = " << mediumMat->GetName() << G4endl;

  // plasma energy squared for plate material
  fSigma3 = fPlasmaCof * mediumMat->GetElectronDensity();
  if(verboseLevel > 0)
    G4cout << "medium plasma energy = " << std::sqrt(fSigma3) / eV << " eV"
           << G4endl;

  // Compute cofs for preparation of linear photo absorption in external medium
  ComputeMediumPhotoAbsCof();
}

///////////////////////////////////////////////////////////////////////////
G4StrawTubeXTRadiator::~G4StrawTubeXTRadiator() = default;

void G4StrawTubeXTRadiator::ProcessDescription(std::ostream& out) const
{
  out << "Simulation of forward X-ray transition radiation for the case of\n"
         "a straw tube radiator.\n";
}

///////////////////////////////////////////////////////////////////////////
// Approximation for radiator interference factor for the case of
// straw tube radiator. The plate (window, straw wall) and gas (inside straw)
// gap thicknesses are gamma distributed.
// The mean values of the plate and gas gap thicknesses
// are supposed to be about XTR formation zone.
G4double G4StrawTubeXTRadiator::GetStackFactor(G4double energy, G4double gamma,
                                               G4double varAngle)
{
  G4double result, L2, L3, M2, M3;

  L2 = GetPlateFormationZone(energy, gamma, varAngle);
  L3 = GetGasFormationZone(energy, gamma, varAngle);

  M2 = GetPlateLinearPhotoAbs(energy);
  M3 = GetGasLinearPhotoAbs(energy);

  G4complex C2(1.0 + 0.5 * fPlateThick * M2 / fAlphaPlate,
               fPlateThick / L2 / fAlphaPlate);
  G4complex C3(1.0 + 0.5 * fGasThick * M3 / fAlphaGas,
               fGasThick / L3 / fAlphaGas);

  G4complex H2 = std::pow(C2, -fAlphaPlate);
  G4complex H3 = std::pow(C3, -fAlphaGas);
  G4complex H  = H2 * H3;

  G4complex Z1 = GetMediumComplexFZ(energy, gamma, varAngle);
  G4complex Z2 = GetPlateComplexFZ(energy, gamma, varAngle);
  G4complex Z3 = GetGasComplexFZ(energy, gamma, varAngle);

  G4complex R = (Z1 - Z2) * (Z1 - Z2) * (1. - H2 * H) +
                (Z2 - Z3) * (Z2 - Z3) * (1. - H3) +
                2. * (Z1 - Z2) * (Z2 - Z3) * H2 * (1. - H3);

  result = 2.0 * std::real(R) * (varAngle * energy / hbarc / hbarc);

  return result;
}

////////////////////////////////////////////////////////////////////////
// Calculates formation zone for external medium. Omega is energy !!!
G4double G4StrawTubeXTRadiator::GetMediumFormationZone(G4double omega,
                                                       G4double gamma,
                                                       G4double varAngle)
{
  G4double cof, lambda;
  lambda = 1.0 / gamma / gamma + varAngle + fSigma3 / omega / omega;
  cof    = 2.0 * hbarc / omega / lambda;
  return cof;
}

////////////////////////////////////////////////////////////////////////
// Calculates complex formation zone for external medium. Omega is energy !!!
G4complex G4StrawTubeXTRadiator::GetMediumComplexFZ(G4double omega,
                                                    G4double gamma,
                                                    G4double varAngle)
{
  G4double cof, length, delta, real_v, image_v;

  length = 0.5 * GetMediumFormationZone(omega, gamma, varAngle);
  delta  = length * GetMediumLinearPhotoAbs(omega);
  cof    = 1.0 / (1.0 + delta * delta);

  real_v  = length * cof;
  image_v = real_v * delta;

  G4complex zone(real_v, image_v);
  return zone;
}

////////////////////////////////////////////////////////////////////////
// Computes matrix of Sandia photo absorption cross section coefficients for
// medium material
void G4StrawTubeXTRadiator::ComputeMediumPhotoAbsCof()
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4Material* mat                   = (*theMaterialTable)[fMatIndex3];
  fMediumPhotoAbsCof                      = mat->GetSandiaTable();
}

//////////////////////////////////////////////////////////////////////
// Returns the value of linear photo absorption coefficient (in reciprocal
// length) for medium for given energy of X-ray photon omega
G4double G4StrawTubeXTRadiator::GetMediumLinearPhotoAbs(G4double omega)
{
  G4double omega2, omega3, omega4;

  omega2 = omega * omega;
  omega3 = omega2 * omega;
  omega4 = omega2 * omega2;

  const G4double* SandiaCof =
    fMediumPhotoAbsCof->GetSandiaCofForMaterial(omega);

  G4double cross = SandiaCof[0] / omega + SandiaCof[1] / omega2 +
                   SandiaCof[2] / omega3 + SandiaCof[3] / omega4;
  return cross;
}
