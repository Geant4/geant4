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
// History:
// 1st version 11.09.97 V. Grichine (Vladimir.Grichine@cern.ch )
// 2nd version 17.12.97 V. Grichine
// 17-09-01, migration of Materials to pure STL (mma)
// 10-03-03, migration to "cut per region" (V.Ivanchenko)
// 03.06.03, V.Ivanchenko fix compilation warnings

#include "G4ForwardXrayTR.hh"

#include "globals.hh"
#include "G4Gamma.hh"
#include "G4GeometryTolerance.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4Poisson.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsModelCatalog.hh"

//////////////////////////////////////////////////////////////////////
//
// Constructor for creation of physics tables (angle and energy TR
// distributions) for a couple of selected materials.
//
// Recommended for use in applications with many materials involved,
// when only few (usually couple) materials are interested for generation
// of TR on the interface between them
G4ForwardXrayTR::G4ForwardXrayTR(const G4String& matName1,
                                 const G4String& matName2,
                                 const G4String& processName)
  : G4TransitionRadiation(processName)
{
  secID = G4PhysicsModelCatalog::GetModelID("model_XrayTR");
  fPtrGamma                = nullptr;
  fGammaCutInKineticEnergy = nullptr;
  fGammaTkinCut = fMinEnergyTR = fMaxEnergyTR = fMaxThetaTR = 0.0;
  fGamma = fSigma1 = fSigma2 = 0.0;
  fAngleDistrTable           = nullptr;
  fEnergyDistrTable          = nullptr;
  fMatIndex1 = fMatIndex2 = 0;

  // Proton energy vector initialization
  fProtonEnergyVector =
    new G4PhysicsLogVector(fMinProtonTkin, fMaxProtonTkin, fTotBin);
  G4int iMat;
  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  G4bool build = true;

  for(iMat = 0; iMat < numOfCouples; ++iMat)  // check first material name
  {
    const G4MaterialCutsCouple* couple =
      theCoupleTable->GetMaterialCutsCouple(iMat);
    if(matName1 == couple->GetMaterial()->GetName())
    {
      fMatIndex1 = couple->GetIndex();
      break;
    }
  }
  if(iMat == numOfCouples)
  {
    G4Exception("G4ForwardXrayTR::G4ForwardXrayTR", "ForwardXrayTR01",
                JustWarning,
                "Invalid first material name in G4ForwardXrayTR constructor!");
    build = false;
  }

  if(build)
  {
    for(iMat = 0; iMat < numOfCouples; ++iMat)  // check second material name
    {
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple(iMat);
      if(matName2 == couple->GetMaterial()->GetName())
      {
        fMatIndex2 = couple->GetIndex();
        break;
      }
    }
    if(iMat == numOfCouples)
    {
      G4Exception(
        "G4ForwardXrayTR::G4ForwardXrayTR", "ForwardXrayTR02", JustWarning,
        "Invalid second material name in G4ForwardXrayTR constructor!");
      build = false;
    }
  }
  if(build)
  {
    BuildXrayTRtables();
  }
}

/////////////////////////////////////////////////////////////////////////
// Constructor used by X-ray transition radiation parametrisation models
G4ForwardXrayTR::G4ForwardXrayTR(const G4String& processName)
  : G4TransitionRadiation(processName)
{
  fPtrGamma                = nullptr;
  fGammaCutInKineticEnergy = nullptr;
  fGammaTkinCut = fMinEnergyTR = fMaxEnergyTR = fMaxThetaTR = 0.0;
  fGamma = fSigma1 = fSigma2 = 0.0;
  fAngleDistrTable           = nullptr;
  fEnergyDistrTable          = nullptr;
  fMatIndex1 = fMatIndex2 = 0;

  // Proton energy vector initialization
  fProtonEnergyVector =
    new G4PhysicsLogVector(fMinProtonTkin, fMaxProtonTkin, fTotBin);
}

//////////////////////////////////////////////////////////////////////
// Destructor
G4ForwardXrayTR::~G4ForwardXrayTR()
{
  delete fAngleDistrTable;
  delete fEnergyDistrTable;
  delete fProtonEnergyVector;
}

void G4ForwardXrayTR::ProcessDescription(std::ostream& out) const
{
  out << "Simulation of forward X-ray transition radiation generated by\n"
         "relativistic charged particles crossing the interface between\n"
         "two materials.\n";
}

G4double G4ForwardXrayTR::GetMeanFreePath(const G4Track&, G4double,
                                          G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;  // so TR doesn't limit mean free path
}

//////////////////////////////////////////////////////////////////////////////
// Build physics tables for energy and angular distributions of X-ray TR photon
void G4ForwardXrayTR::BuildXrayTRtables()
{
  G4int iMat, jMat, iTkin, iTR, iPlace;
  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  fGammaCutInKineticEnergy = theCoupleTable->GetEnergyCutsVector(idxG4GammaCut);

  fAngleDistrTable  = new G4PhysicsTable(2 * fTotBin);
  fEnergyDistrTable = new G4PhysicsTable(2 * fTotBin);

  for(iMat = 0; iMat < numOfCouples;
      ++iMat)  // loop over pairs of different materials
  {
    if(iMat != fMatIndex1 && iMat != fMatIndex2)
      continue;

    for(jMat = 0; jMat < numOfCouples; ++jMat)  // transition iMat -> jMat !!!
    {
      if(iMat == jMat || (jMat != fMatIndex1 && jMat != fMatIndex2))
      {
        continue;
      }
      else
      {
        const G4MaterialCutsCouple* iCouple =
          theCoupleTable->GetMaterialCutsCouple(iMat);
        const G4MaterialCutsCouple* jCouple =
          theCoupleTable->GetMaterialCutsCouple(jMat);
        const G4Material* mat1 = iCouple->GetMaterial();
        const G4Material* mat2 = jCouple->GetMaterial();

        fSigma1 = fPlasmaCof * (mat1->GetElectronDensity());
        fSigma2 = fPlasmaCof * (mat2->GetElectronDensity());

        fGammaTkinCut = 0.0;

        if(fGammaTkinCut > fTheMinEnergyTR)  // setting of min/max TR energies
        {
          fMinEnergyTR = fGammaTkinCut;
        }
        else
        {
          fMinEnergyTR = fTheMinEnergyTR;
        }
        if(fGammaTkinCut > fTheMaxEnergyTR)
        {
          fMaxEnergyTR = 2.0 * fGammaTkinCut;  // usually very low TR rate
        }
        else
        {
          fMaxEnergyTR = fTheMaxEnergyTR;
        }
        for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
        {
          auto energyVector =
            new G4PhysicsLogVector(fMinEnergyTR, fMaxEnergyTR, fBinTR);

          fGamma = 1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) /
                          proton_mass_c2);

          fMaxThetaTR = 10000.0 / (fGamma * fGamma);

          if(fMaxThetaTR > fTheMaxAngle)
          {
            fMaxThetaTR = fTheMaxAngle;
          }
          else
          {
            if(fMaxThetaTR < fTheMinAngle)
            {
              fMaxThetaTR = fTheMinAngle;
            }
          }
          auto angleVector =
            new G4PhysicsLinearVector(0.0, fMaxThetaTR, fBinTR);
          G4double energySum = 0.0;
          G4double angleSum  = 0.0;

          energyVector->PutValue(fBinTR - 1, energySum);
          angleVector->PutValue(fBinTR - 1, angleSum);

          for(iTR = fBinTR - 2; iTR >= 0; --iTR)
          {
            energySum +=
              fCofTR * EnergySum(energyVector->GetLowEdgeEnergy(iTR),
                                 energyVector->GetLowEdgeEnergy(iTR + 1));

            angleSum +=
              fCofTR * AngleSum(angleVector->GetLowEdgeEnergy(iTR),
                                angleVector->GetLowEdgeEnergy(iTR + 1));

            energyVector->PutValue(iTR, energySum);
            angleVector->PutValue(iTR, angleSum);
          }

          if(jMat < iMat)
          {
            iPlace = fTotBin + iTkin;
          }
          else  // jMat > iMat right part of matrices (jMat-1) !
          {
            iPlace = iTkin;
          }
          fEnergyDistrTable->insertAt(iPlace, energyVector);
          fAngleDistrTable->insertAt(iPlace, angleVector);
        }  //                      iTkin
      }    //         jMat != iMat
    }      //     jMat
  }        // iMat
}

///////////////////////////////////////////////////////////////////////
//
// This function returns the spectral and angle density of TR quanta
// in X-ray energy region generated forward when a relativistic
// charged particle crosses interface between two materials.
// The high energy small theta approximation is applied.
// (matter1 -> matter2)
// varAngle =2* (1 - std::cos(Theta)) or approximately = Theta*Theta
//
G4double G4ForwardXrayTR::SpectralAngleTRdensity(G4double energy,
                                                 G4double varAngle) const
{
  G4double formationLength1, formationLength2;
  formationLength1 =
    1.0 / (1.0 / (fGamma * fGamma) + fSigma1 / (energy * energy) + varAngle);
  formationLength2 =
    1.0 / (1.0 / (fGamma * fGamma) + fSigma2 / (energy * energy) + varAngle);
  return (varAngle / energy) * (formationLength1 - formationLength2) *
         (formationLength1 - formationLength2);
}

//////////////////////////////////////////////////////////////////
// Analytical formula for angular density of X-ray TR photons
G4double G4ForwardXrayTR::AngleDensity(G4double energy, G4double varAngle) const
{
  G4double x, x2, c, d, f, a2, b2, a4, b4;
  G4double cof1, cof2, cof3;
  x    = 1.0 / energy;
  x2   = x * x;
  c    = 1.0 / fSigma1;
  d    = 1.0 / fSigma2;
  f    = (varAngle + 1.0 / (fGamma * fGamma));
  a2   = c * f;
  b2   = d * f;
  a4   = a2 * a2;
  b4   = b2 * b2;
  cof1 = c * c * (0.5 / (a2 * (x2 + a2)) + 0.5 * std::log(x2 / (x2 + a2)) / a4);
  cof3 = d * d * (0.5 / (b2 * (x2 + b2)) + 0.5 * std::log(x2 / (x2 + b2)) / b4);
  cof2 = -c * d *
         (std::log(x2 / (x2 + b2)) / b2 - std::log(x2 / (x2 + a2)) / a2) /
         (a2 - b2);
  return -varAngle * (cof1 + cof2 + cof3);
}

/////////////////////////////////////////////////////////////////////
// Definite integral of X-ray TR spectral-angle density from energy1
// to energy2
G4double G4ForwardXrayTR::EnergyInterval(G4double energy1, G4double energy2,
                                         G4double varAngle) const
{
  return AngleDensity(energy2, varAngle) - AngleDensity(energy1, varAngle);
}

//////////////////////////////////////////////////////////////////////
// Integral angle distribution of X-ray TR photons based on analytical
// formula for angle density
G4double G4ForwardXrayTR::AngleSum(G4double varAngle1, G4double varAngle2) const
{
  G4int i;
  G4double h, sumEven = 0.0, sumOdd = 0.0;
  h = 0.5 * (varAngle2 - varAngle1) / fSympsonNumber;
  for(i = 1; i < fSympsonNumber; ++i)
  {
    sumEven +=
      EnergyInterval(fMinEnergyTR, fMaxEnergyTR, varAngle1 + 2 * i * h);
    sumOdd +=
      EnergyInterval(fMinEnergyTR, fMaxEnergyTR, varAngle1 + (2 * i - 1) * h);
  }
  sumOdd += EnergyInterval(fMinEnergyTR, fMaxEnergyTR,
                           varAngle1 + (2 * fSympsonNumber - 1) * h);

  return h *
         (EnergyInterval(fMinEnergyTR, fMaxEnergyTR, varAngle1) +
          EnergyInterval(fMinEnergyTR, fMaxEnergyTR, varAngle2) + 4.0 * sumOdd +
          2.0 * sumEven) /
         3.0;
}

/////////////////////////////////////////////////////////////////////
// Analytical Expression for   spectral density of Xray TR photons
// x = 2*(1 - std::cos(Theta)) ~ Theta^2
G4double G4ForwardXrayTR::SpectralDensity(G4double energy, G4double x) const
{
  G4double a, b;
  a = 1.0 / (fGamma * fGamma) + fSigma1 / (energy * energy);
  b = 1.0 / (fGamma * fGamma) + fSigma2 / (energy * energy);
  return ((a + b) * std::log((x + b) / (x + a)) / (a - b) + a / (x + a) +
          b / (x + b)) /
         energy;
}

////////////////////////////////////////////////////////////////////
//  The spectral density in some angle interval from varAngle1 to
//  varAngle2
G4double G4ForwardXrayTR::AngleInterval(G4double energy, G4double varAngle1,
                                        G4double varAngle2) const
{
  return SpectralDensity(energy, varAngle2) -
         SpectralDensity(energy, varAngle1);
}

////////////////////////////////////////////////////////////////////
// Integral spectral distribution of X-ray TR photons based on
// analytical formula for spectral density
G4double G4ForwardXrayTR::EnergySum(G4double energy1, G4double energy2) const
{
  G4int i;
  G4double h, sumEven = 0.0, sumOdd = 0.0;
  h = 0.5 * (energy2 - energy1) / fSympsonNumber;
  for(i = 1; i < fSympsonNumber; ++i)
  {
    sumEven += AngleInterval(energy1 + 2 * i * h, 0.0, fMaxThetaTR);
    sumOdd += AngleInterval(energy1 + (2 * i - 1) * h, 0.0, fMaxThetaTR);
  }
  sumOdd +=
    AngleInterval(energy1 + (2 * fSympsonNumber - 1) * h, 0.0, fMaxThetaTR);

  return h *
         (AngleInterval(energy1, 0.0, fMaxThetaTR) +
          AngleInterval(energy2, 0.0, fMaxThetaTR) + 4.0 * sumOdd +
          2.0 * sumEven) /
         3.0;
}

/////////////////////////////////////////////////////////////////////////
// PostStepDoIt function for creation of forward X-ray photons in TR process
// on boundary between two materials with really different plasma energies
G4VParticleChange* G4ForwardXrayTR::PostStepDoIt(const G4Track& aTrack,
                                                 const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  G4int iMat, jMat, iTkin, iPlace, numOfTR, iTR, iTransfer;

  G4double energyPos, anglePos, energyTR, theta, phi, dirX, dirY, dirZ;
  G4double W, W1, W2, E1, E2;

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4double tol =
    0.5 * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if(pPostStepPoint->GetStepStatus() != fGeomBoundary)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  if(aTrack.GetStepLength() <= tol)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  // Arrived at boundary, so begin to try TR

  const G4MaterialCutsCouple* iCouple = pPreStepPoint->GetPhysicalVolume()
                                          ->GetLogicalVolume()
                                          ->GetMaterialCutsCouple();
  const G4MaterialCutsCouple* jCouple = pPostStepPoint->GetPhysicalVolume()
                                          ->GetLogicalVolume()
                                          ->GetMaterialCutsCouple();
  const G4Material* iMaterial = iCouple->GetMaterial();
  const G4Material* jMaterial = jCouple->GetMaterial();
  iMat                        = iCouple->GetIndex();
  jMat                        = jCouple->GetIndex();

  // The case of equal or approximate (in terms of plasma energy) materials
  // No TR photons ?!

  if(iMat == jMat ||
     ((fMatIndex1 >= 0 && fMatIndex2 >= 0) &&
      (iMat != fMatIndex1 && iMat != fMatIndex2) &&
      (jMat != fMatIndex1 && jMat != fMatIndex2))

     || iMaterial->GetState() == jMaterial->GetState()

     || (iMaterial->GetState() == kStateSolid &&
         jMaterial->GetState() == kStateLiquid)

     || (iMaterial->GetState() == kStateLiquid &&
         jMaterial->GetState() == kStateSolid))
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double charge = aParticle->GetDefinition()->GetPDGCharge();

  if(charge == 0.0)  // Uncharged particle doesn't Generate TR photons
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  // Now we are ready to Generate TR photons

  G4double chargeSq  = charge * charge;
  G4double kinEnergy = aParticle->GetKineticEnergy();
  G4double massRatio =
    proton_mass_c2 / aParticle->GetDefinition()->GetPDGMass();
  G4double TkinScaled = kinEnergy * massRatio;
  for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  {
    if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))
    {
      break;
    }
  }
  if(jMat < iMat)
  {
    iPlace = fTotBin + iTkin - 1;
  }
  else
  {
    iPlace = iTkin - 1;
  }

  G4ParticleMomentum particleDir = aParticle->GetMomentumDirection();

  if(iTkin == fTotBin)  // TR plato, try from left
  {
    numOfTR = (G4int)G4Poisson(
      ((*(*fEnergyDistrTable)(iPlace))(0) + (*(*fAngleDistrTable)(iPlace))(0)) *
      chargeSq * 0.5);
    if(numOfTR == 0)
    {
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    else
    {
      aParticleChange.SetNumberOfSecondaries(numOfTR);

      for(iTR = 0; iTR < numOfTR; ++iTR)
      {
        energyPos = (*(*fEnergyDistrTable)(iPlace))(0) * G4UniformRand();
        for(iTransfer = 0; iTransfer < fBinTR - 1; ++iTransfer)
        {
          if(energyPos >= (*(*fEnergyDistrTable)(iPlace))(iTransfer))
            break;
        }
        energyTR = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

        kinEnergy -= energyTR;
        aParticleChange.ProposeEnergy(kinEnergy);

        anglePos = (*(*fAngleDistrTable)(iPlace))(0) * G4UniformRand();
        for(iTransfer = 0; iTransfer < fBinTR - 1; ++iTransfer)
        {
          if(anglePos > (*(*fAngleDistrTable)(iPlace))(iTransfer))
            break;
        }
        theta = std::sqrt(
          (*fAngleDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer - 1));

        phi  = twopi * G4UniformRand();
        dirX = std::sin(theta) * std::cos(phi);
        dirY = std::sin(theta) * std::sin(phi);
        dirZ = std::cos(theta);
        G4ThreeVector directionTR(dirX, dirY, dirZ);
        directionTR.rotateUz(particleDir);
        auto aPhotonTR = new G4DynamicParticle(G4Gamma::Gamma(), directionTR, energyTR);

	// Create the G4Track
	auto aSecondaryTrack = new G4Track(aPhotonTR, aTrack.GetGlobalTime(), aTrack.GetPosition());
	aSecondaryTrack->SetTouchableHandle(aStep.GetPostStepPoint()->GetTouchableHandle());
	aSecondaryTrack->SetParentID(aTrack.GetTrackID());
	aSecondaryTrack->SetCreatorModelID(secID);
	aParticleChange.AddSecondary(aSecondaryTrack);	
      }
    }
  }
  else
  {
    if(iTkin == 0)  // Tkin is too small, neglect of TR photon generation
    {
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    else  // general case: Tkin between two vectors of the material
    {
      E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1);
      E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
      W  = 1.0 / (E2 - E1);
      W1 = (E2 - TkinScaled) * W;
      W2 = (TkinScaled - E1) * W;

      numOfTR = (G4int)G4Poisson((((*(*fEnergyDistrTable)(iPlace))(0) +
                                   (*(*fAngleDistrTable)(iPlace))(0)) *
                                    W1 +
                                  ((*(*fEnergyDistrTable)(iPlace + 1))(0) +
                                   (*(*fAngleDistrTable)(iPlace + 1))(0)) *
                                    W2) *
                                 chargeSq * 0.5);
      if(numOfTR == 0)
      {
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
      else
      {
        aParticleChange.SetNumberOfSecondaries(numOfTR);
        for(iTR = 0; iTR < numOfTR; ++iTR)
        {
          energyPos = ((*(*fEnergyDistrTable)(iPlace))(0) * W1 +
                       (*(*fEnergyDistrTable)(iPlace + 1))(0) * W2) *
                      G4UniformRand();
          for(iTransfer = 0; iTransfer < fBinTR - 1; ++iTransfer)
          {
            if(energyPos >=
               ((*(*fEnergyDistrTable)(iPlace))(iTransfer) *W1 +
                (*(*fEnergyDistrTable)(iPlace + 1))(iTransfer) *W2))
              break;
          }
          energyTR =
            ((*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer)) * W1 +
            ((*fEnergyDistrTable)(iPlace + 1)->GetLowEdgeEnergy(iTransfer)) *
              W2;

          kinEnergy -= energyTR;
          aParticleChange.ProposeEnergy(kinEnergy);

          anglePos = ((*(*fAngleDistrTable)(iPlace))(0) * W1 +
                      (*(*fAngleDistrTable)(iPlace + 1))(0) * W2) *
                     G4UniformRand();
          for(iTransfer = 0; iTransfer < fBinTR - 1; ++iTransfer)
          {
            if(anglePos > ((*(*fAngleDistrTable)(iPlace))(iTransfer) *W1 +
                           (*(*fAngleDistrTable)(iPlace + 1))(iTransfer) *W2))
              break;
          }
          theta = std::sqrt(
            ((*fAngleDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer - 1)) *
              W1 +
            ((*fAngleDistrTable)(iPlace + 1)->GetLowEdgeEnergy(iTransfer - 1)) *
              W2);

          phi  = twopi * G4UniformRand();
          dirX = std::sin(theta) * std::cos(phi);
          dirY = std::sin(theta) * std::sin(phi);
          dirZ = std::cos(theta);
          G4ThreeVector directionTR(dirX, dirY, dirZ);
          directionTR.rotateUz(particleDir);
          auto aPhotonTR =
            new G4DynamicParticle(G4Gamma::Gamma(), directionTR, energyTR);

	  // Create the G4Track
	  G4Track* aSecondaryTrack = new G4Track(aPhotonTR, aTrack.GetGlobalTime(), aTrack.GetPosition());
	  aSecondaryTrack->SetTouchableHandle(aStep.GetPostStepPoint()->GetTouchableHandle());
	  aSecondaryTrack->SetParentID(aTrack.GetTrackID());
	  aSecondaryTrack->SetCreatorModelID(secID);
	  aParticleChange.AddSecondary(aSecondaryTrack);		  
        }
      }
    }
  }
  return &aParticleChange;
}

////////////////////////////////////////////////////////////////////////////
// Test function for checking of PostStepDoIt random preparation of TR photon
// energy
G4double G4ForwardXrayTR::GetEnergyTR(G4int iMat, G4int jMat, G4int iTkin) const
{
  G4int iPlace, numOfTR, iTR, iTransfer;
  G4double energyTR = 0.0;  // return this value for no TR photons
  G4double energyPos;
  G4double W1, W2;

  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  // The case of equal or approximate (in terms of plasma energy) materials
  // No TR photons ?!

  const G4MaterialCutsCouple* iCouple =
    theCoupleTable->GetMaterialCutsCouple(iMat);
  const G4MaterialCutsCouple* jCouple =
    theCoupleTable->GetMaterialCutsCouple(jMat);
  const G4Material* iMaterial = iCouple->GetMaterial();
  const G4Material* jMaterial = jCouple->GetMaterial();

  if(iMat == jMat

     || iMaterial->GetState() == jMaterial->GetState()

     || (iMaterial->GetState() == kStateSolid &&
         jMaterial->GetState() == kStateLiquid)

     || (iMaterial->GetState() == kStateLiquid &&
         jMaterial->GetState() == kStateSolid))

  {
    return energyTR;
  }

  if(jMat < iMat)
  {
    iPlace = (iMat * (numOfCouples - 1) + jMat) * fTotBin + iTkin - 1;
  }
  else
  {
    iPlace = (iMat * (numOfCouples - 1) + jMat - 1) * fTotBin + iTkin - 1;
  }
  G4PhysicsVector* energyVector1 = (*fEnergyDistrTable)(iPlace);
  G4PhysicsVector* energyVector2 = (*fEnergyDistrTable)(iPlace + 1);

  if(iTkin == fTotBin)  // TR plato, try from left
  {
    numOfTR = (G4int)G4Poisson((*energyVector1)(0));
    if(numOfTR == 0)
    {
      return energyTR;
    }
    else
    {
      for(iTR = 0; iTR < numOfTR; ++iTR)
      {
        energyPos = (*energyVector1)(0) * G4UniformRand();
        for(iTransfer = 0; iTransfer < fBinTR - 1; ++iTransfer)
        {
          if(energyPos >= (*energyVector1)(iTransfer))
            break;
        }
        energyTR += energyVector1->GetLowEdgeEnergy(iTransfer);
      }
    }
  }
  else
  {
    if(iTkin == 0)  // Tkin is too small, neglect of TR photon generation
    {
      return energyTR;
    }
    else  // general case: Tkin between two vectors of the material
    {     // use trivial mean half/half
      W1      = 0.5;
      W2      = 0.5;
      numOfTR = (G4int)G4Poisson((*energyVector1)(0) * W1 + (*energyVector2)(0) * W2);
      if(numOfTR == 0)
      {
        return energyTR;
      }
      else
      {
        G4cout << "It is still OK in GetEnergyTR(int,int,int)" << G4endl;
        for(iTR = 0; iTR < numOfTR; ++iTR)
        {
          energyPos = ((*energyVector1)(0) * W1 + (*energyVector2)(0) * W2) *
                      G4UniformRand();
          for(iTransfer = 0; iTransfer < fBinTR - 1; ++iTransfer)
          {
            if(energyPos >= ((*energyVector1)(iTransfer) *W1 +
                             (*energyVector2)(iTransfer) *W2))
              break;
          }
          energyTR += (energyVector1->GetLowEdgeEnergy(iTransfer)) * W1 +
                      (energyVector2->GetLowEdgeEnergy(iTransfer)) * W2;
        }
      }
    }
  }

  return energyTR;
}

////////////////////////////////////////////////////////////////////////////
// Test function for checking of PostStepDoIt random preparation of TR photon
// theta angle relative to particle direction
G4double G4ForwardXrayTR::GetThetaTR(G4int, G4int, G4int) const { return 0.0; }

G4int G4ForwardXrayTR::GetSympsonNumber() { return fSympsonNumber; }

G4int G4ForwardXrayTR::GetBinTR() { return fBinTR; }

G4double G4ForwardXrayTR::GetMinProtonTkin() { return fMinProtonTkin; }

G4double G4ForwardXrayTR::GetMaxProtonTkin() { return fMaxProtonTkin; }

G4int G4ForwardXrayTR::GetTotBin() { return fTotBin; }
