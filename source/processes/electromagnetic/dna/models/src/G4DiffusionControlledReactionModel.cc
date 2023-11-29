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
//
// Author: Hoang TRAN

#include "G4DiffusionControlledReactionModel.hh"
#include "G4Track.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Exp.hh"
#include "G4IRTUtils.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron_aq.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"
#include "G4ErrorFunction.hh"

G4DiffusionControlledReactionModel::G4DiffusionControlledReactionModel()
  : G4VDNAReactionModel()
{}

G4DiffusionControlledReactionModel::~G4DiffusionControlledReactionModel() =
  default;

void G4DiffusionControlledReactionModel::Initialise(
  const G4MolecularConfiguration* pMolecule, const G4Track&)
{
  fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

void G4DiffusionControlledReactionModel::InitialiseToPrint(
  const G4MolecularConfiguration* pMolecule)
{
  fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

G4double G4DiffusionControlledReactionModel::GetReactionRadius(
  const G4MolecularConfiguration* pMol1, const G4MolecularConfiguration* pMol2)
{
  auto reactionData = fpReactionTable->GetReactionData(pMol1, pMol2);
  if(reactionData == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No reactionData"
                         << " for : " << pMol1->GetName() << " and "
                         << pMol2->GetName();
    G4Exception("G4DiffusionControlledReactionModel"
                "::GetReactionRadius()",
                "G4DiffusionControlledReactionModel00", FatalException,
                exceptionDescription);
    return 0.;
  }else
  {
    return reactionData->GetEffectiveReactionRadius();
  }
}

G4double G4DiffusionControlledReactionModel::GetReactionRadius(const G4int& i)
{
  auto pMol1 = (*fpReactionData)[i]->GetReactant1();
  auto pMol2 = (*fpReactionData)[i]->GetReactant2();
  return GetReactionRadius(pMol1, pMol2);
}

G4double G4DiffusionControlledReactionModel::GetTimeToEncounter(
  const G4Track& trackA, const G4Track& trackB)
{
  auto pMolConfA = GetMolecule(trackA)->GetMolecularConfiguration();
  auto pMolConfB = GetMolecule(trackB)->GetMolecularConfiguration();

  G4double D =
    pMolConfA->GetDiffusionCoefficient() + pMolConfB->GetDiffusionCoefficient();

  if(D == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "The total diffusion coefficient for : "
                         << pMolConfA->GetName() << " and "
                         << pMolConfB->GetName() << " is null ";
    G4Exception("G4DiffusionControlledReactionModel"
                "::GetTimeToEncounter()",
                "G4DiffusionControlledReactionModel03", FatalException,
                exceptionDescription);
  }

  auto reactionData = G4DNAMolecularReactionTable::Instance()->GetReactionData(
    pMolConfA, pMolConfB);
  G4double kobs     = reactionData->GetObservedReactionRateConstant();
  G4double distance = (trackA.GetPosition() - trackB.GetPosition()).mag();
  G4double SmoluchowskiRadius = reactionData->GetEffectiveReactionRadius();

  if(distance == 0 || distance < SmoluchowskiRadius)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "distance = " << distance << " is uncorrected with "
                         << " Reff = " << SmoluchowskiRadius
                         << " for : " << pMolConfA->GetName() << " and "
                         << pMolConfB->GetName();
    G4Exception("G4DiffusionControlledReactionModel"
                "::GetTimeToEncounter()",
                "G4DiffusionControlledReactionModel02", FatalException,
                exceptionDescription);
  }
  else
  {
    G4double Winf  = SmoluchowskiRadius / distance;
    G4double U     = G4UniformRand();
    G4double X     = 0;
    G4double irt_1 = -1.0 * ps;
    if(Winf > 0 && U < Winf)
    {
      G4double erfcIn = G4ErrorFunction::erfcInv(U / Winf);
      if(erfcIn != 0)
      {
        G4double d =
          (distance - SmoluchowskiRadius) / erfcIn;
        irt_1 = (1.0 / (4 * D)) * d * d;
      }
    }

    if(reactionData->GetReactionType() == 0)  // Totally diffused contr
    {
      return irt_1;
    }

    if(irt_1 < 0)
    {
      return irt_1;
    }
    else
    {
      G4double kdif = 4 * CLHEP::pi * D * SmoluchowskiRadius * Avogadro;

      if(pMolConfA == pMolConfB)
      {
        kdif /= 2;
      }
      G4double kact   = G4IRTUtils::GetKact(kobs, kdif);
      G4double sumOfk = kact + kdif;
      if(sumOfk != 0)
      {
        G4double rateFactor = kact / sumOfk;
        if(G4UniformRand() > rateFactor)
        {
          return -1.0 * ps;
        }
        G4double Y = std::abs(G4RandGauss::shoot(0.0, std::sqrt(2)));

        if(Y > 0)
        {
          X = -(G4Log(G4UniformRand())) / Y;
        }
        G4double f     = X * SmoluchowskiRadius * kdif / sumOfk;
        G4double irt_2 = (f * f) / D;
        return irt_1 + irt_2;
      }
    }
  }
  return -1.0 * ps;
}
