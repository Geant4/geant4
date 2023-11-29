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
// G4TransitionRadiation class -- implementation file

// GEANT 4 class implementation file --- Copyright CERN 1995

// For information related to this code, please, contact
// CERN, CN Division, ASD Group
// History:
// 1st version 11.09.97 V. Grichine (Vladimir.Grichine@cern.ch )
// 2nd version 16.12.97 V. Grichine
// 3rd version 28.07.05, P.Gumplinger add G4ProcessType to constructor

//#include <cmath>

#include "G4TransitionRadiation.hh"

#include "G4EmProcessSubType.hh"

///////////////////////////////////////////////////////////////////////
// Constructor for selected couple of materials
G4TransitionRadiation::G4TransitionRadiation(const G4String& processName,
                                             G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  SetProcessSubType(fTransitionRadiation);
  fMatIndex1 = fMatIndex2 = 0;

  fGamma = fEnergy = fVarAngle = fMinEnergy = fMaxEnergy = fMaxTheta = 0.0;
  fSigma1 = fSigma2 = 0.0;
}

//////////////////////////////////////////////////////////////////////
// Destructor
G4TransitionRadiation::~G4TransitionRadiation() = default;

void G4TransitionRadiation::ProcessDescription(std::ostream& out) const
{
  out << "Base class for simulation of x-ray transition radiation.\n";
}

G4bool G4TransitionRadiation::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return (aParticleType.GetPDGCharge() != 0.0);
}

G4double G4TransitionRadiation::GetMeanFreePath(const G4Track&, G4double,
                                                G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;  // so TR doesn't limit mean free path
}

G4VParticleChange* G4TransitionRadiation::PostStepDoIt(const G4Track&,
                                                       const G4Step&)
{
  ClearNumberOfInteractionLengthLeft();
  return &aParticleChange;
}

///////////////////////////////////////////////////////////////////
// Sympson integral of TR spectral-angle density over energy between
// the limits energy 1 and energy2 at fixed varAngle = 1 - std::cos(Theta)
G4double G4TransitionRadiation::IntegralOverEnergy(G4double energy1,
                                                   G4double energy2,
                                                   G4double varAngle) const
{
  G4int i;
  G4double h, sumEven = 0.0, sumOdd = 0.0;
  h = 0.5 * (energy2 - energy1) / fSympsonNumber;
  for(i = 1; i < fSympsonNumber; i++)
  {
    sumEven += SpectralAngleTRdensity(energy1 + 2 * i * h, varAngle);
    sumOdd += SpectralAngleTRdensity(energy1 + (2 * i - 1) * h, varAngle);
  }
  sumOdd +=
    SpectralAngleTRdensity(energy1 + (2 * fSympsonNumber - 1) * h, varAngle);
  return h *
         (SpectralAngleTRdensity(energy1, varAngle) +
          SpectralAngleTRdensity(energy2, varAngle) + 4.0 * sumOdd +
          2.0 * sumEven) /
         3.0;
}

///////////////////////////////////////////////////////////////////
// Sympson integral of TR spectral-angle density over energy between
// the limits varAngle1 and varAngle2 at fixed energy
G4double G4TransitionRadiation::IntegralOverAngle(G4double energy,
                                                  G4double varAngle1,
                                                  G4double varAngle2) const
{
  G4int i;
  G4double h, sumEven = 0.0, sumOdd = 0.0;
  h = 0.5 * (varAngle2 - varAngle1) / fSympsonNumber;
  for(i = 1; i < fSympsonNumber; ++i)
  {
    sumEven += SpectralAngleTRdensity(energy, varAngle1 + 2 * i * h);
    sumOdd += SpectralAngleTRdensity(energy, varAngle1 + (2 * i - 1) * h);
  }
  sumOdd +=
    SpectralAngleTRdensity(energy, varAngle1 + (2 * fSympsonNumber - 1) * h);

  return h *
         (SpectralAngleTRdensity(energy, varAngle1) +
          SpectralAngleTRdensity(energy, varAngle2) + 4.0 * sumOdd +
          2.0 * sumEven) /
         3.0;
}

///////////////////////////////////////////////////////////////////
// The number of transition radiation photons generated in the
// angle interval between varAngle1 and varAngle2
G4double G4TransitionRadiation::AngleIntegralDistribution(
  G4double varAngle1, G4double varAngle2) const
{
  G4int i;
  G4double h, sumEven = 0.0, sumOdd = 0.0;
  h = 0.5 * (varAngle2 - varAngle1) / fSympsonNumber;
  for(i = 1; i < fSympsonNumber; ++i)
  {
    sumEven += IntegralOverEnergy(fMinEnergy,
                                  fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                                  varAngle1 + 2 * i * h) +
               IntegralOverEnergy(fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                                  fMaxEnergy, varAngle1 + 2 * i * h);
    sumOdd += IntegralOverEnergy(fMinEnergy,
                                 fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                                 varAngle1 + (2 * i - 1) * h) +
              IntegralOverEnergy(fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                                 fMaxEnergy, varAngle1 + (2 * i - 1) * h);
  }
  sumOdd +=
    IntegralOverEnergy(fMinEnergy, fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                       varAngle1 + (2 * fSympsonNumber - 1) * h) +
    IntegralOverEnergy(fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy), fMaxEnergy,
                       varAngle1 + (2 * fSympsonNumber - 1) * h);

  return h *
         (IntegralOverEnergy(fMinEnergy,
                             fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                             varAngle1) +
          IntegralOverEnergy(fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                             fMaxEnergy, varAngle1) +
          IntegralOverEnergy(fMinEnergy,
                             fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                             varAngle2) +
          IntegralOverEnergy(fMinEnergy + 0.3 * (fMaxEnergy - fMinEnergy),
                             fMaxEnergy, varAngle2) +
          4.0 * sumOdd + 2.0 * sumEven) /
         3.0;
}

///////////////////////////////////////////////////////////////////
// The number of transition radiation photons, generated in the
// energy interval between energy1 and energy2
G4double G4TransitionRadiation::EnergyIntegralDistribution(
  G4double energy1, G4double energy2) const
{
  G4int i;
  G4double h, sumEven = 0.0, sumOdd = 0.0;
  h = 0.5 * (energy2 - energy1) / fSympsonNumber;
  for(i = 1; i < fSympsonNumber; ++i)
  {
    sumEven +=
      IntegralOverAngle(energy1 + 2 * i * h, 0.0, 0.01 * fMaxTheta) +
      IntegralOverAngle(energy1 + 2 * i * h, 0.01 * fMaxTheta, fMaxTheta);
    sumOdd +=
      IntegralOverAngle(energy1 + (2 * i - 1) * h, 0.0, 0.01 * fMaxTheta) +
      IntegralOverAngle(energy1 + (2 * i - 1) * h, 0.01 * fMaxTheta, fMaxTheta);
  }
  sumOdd += IntegralOverAngle(energy1 + (2 * fSympsonNumber - 1) * h, 0.0,
                              0.01 * fMaxTheta) +
            IntegralOverAngle(energy1 + (2 * fSympsonNumber - 1) * h,
                              0.01 * fMaxTheta, fMaxTheta);

  return h *
         (IntegralOverAngle(energy1, 0.0, 0.01 * fMaxTheta) +
          IntegralOverAngle(energy1, 0.01 * fMaxTheta, fMaxTheta) +
          IntegralOverAngle(energy2, 0.0, 0.01 * fMaxTheta) +
          IntegralOverAngle(energy2, 0.01 * fMaxTheta, fMaxTheta) +
          4.0 * sumOdd + 2.0 * sumEven) /
         3.0;
}
