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
// $Id: G4DNASmoluchowskiReactionModel.cc 94218 2015-11-09 08:24:48Z gcosmo $
//
#include "G4DNASmoluchowskiReactionModel.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4UnitsTable.hh"
#include "G4Molecule.hh"
//#include "G4Scheduler.hh"

G4DNASmoluchowskiReactionModel::G4DNASmoluchowskiReactionModel() :
    G4VDNAReactionModel()
{
  fReactionData = 0;
}

G4DNASmoluchowskiReactionModel::G4DNASmoluchowskiReactionModel(const G4DNASmoluchowskiReactionModel& __right) :
    G4VDNAReactionModel(__right)
{
  fReactionData = 0;
}

G4DNASmoluchowskiReactionModel& G4DNASmoluchowskiReactionModel::operator=(const G4DNASmoluchowskiReactionModel& right)
{
  if (this == &right) return *this;
  fReactionData = 0;
  return *this;
}

G4DNASmoluchowskiReactionModel::~G4DNASmoluchowskiReactionModel()
{
  fReactionData = 0;
}

void G4DNASmoluchowskiReactionModel::Initialise(G4MolecularConfiguration* __molecule,
                                                const G4Track&)
{
  fReactionData = fReactionTable->GetReactionData(__molecule);
}

void
G4DNASmoluchowskiReactionModel::
InitialiseToPrint(G4MolecularConfiguration* __molecule)
{
  fReactionData = fReactionTable->GetReactionData(__molecule);
}

G4double
G4DNASmoluchowskiReactionModel::GetReactionRadius(G4MolecularConfiguration* __mol1,
                                                  G4MolecularConfiguration* __mol2)
{
  G4double __output = fReactionTable->GetReactionData(__mol1, __mol2)
      ->GetEffectiveReactionRadius();
  return __output;
}

G4double G4DNASmoluchowskiReactionModel::GetReactionRadius(const G4int __i)
{
  G4double __output = (*fReactionData)[__i]->GetEffectiveReactionRadius();
  return __output;
}

G4bool G4DNASmoluchowskiReactionModel::FindReaction(const G4Track& __trackA,
                                                    const G4Track& __trackB,
                                                    const G4double __R,
                                                    G4double& __r,
                                                    const G4bool __alongStepReaction)
{
  G4double postStepSeparation = 0;
  bool do_break = false;
  G4double R2 = __R * __R;
  int k = 0;

  for (; k < 3; k++)
  {
    postStepSeparation += std::pow(
        __trackA.GetPosition()[k] - __trackB.GetPosition()[k], 2);

    if (postStepSeparation > R2)
    {
      do_break = true;
      break;
    }
  }

  if (do_break == false)
  {
    // The loop was not break
    // => __r^2 < __R^2
    __r = std::sqrt(postStepSeparation);
    return true;
  }
  else if (__alongStepReaction == true)
  {
    //G4cout << "alongStepReaction==true" << G4endl;
    //Along step cheack and
    // the loop has break

    // Continue loop
    for (; k < 3; k++)
    {
      postStepSeparation += std::pow(
          __trackA.GetPosition()[k] - __trackB.GetPosition()[k], 2);
    }
    // Use Green approach : the Brownian bridge
    __r = (postStepSeparation = std::sqrt(postStepSeparation));

    G4Molecule* __moleculeA = GetMolecule(__trackA);
    G4Molecule* __moleculeB = GetMolecule(__trackB);

    G4double __D = __moleculeA->GetDiffusionCoefficient()
        + __moleculeB->GetDiffusionCoefficient();

    G4ThreeVector __preStepPositionA = __trackA.GetStep()->GetPreStepPoint()
        ->GetPosition();
    G4ThreeVector __preStepPositionB = __trackB.GetStep()->GetPreStepPoint()
        ->GetPosition();

    if (__preStepPositionA == __trackA.GetPosition())
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "The molecule : " << __moleculeA->GetName();
      exceptionDescription << " with track ID :" << __trackA.GetTrackID();
      exceptionDescription << " did not move since the previous step." << G4endl;
      exceptionDescription << "Current position : "
                           << G4BestUnit(__trackA.GetPosition(), "Length")
                           << G4endl;
      exceptionDescription << "Previous position : "
                           << G4BestUnit(__preStepPositionA, "Length") << G4endl;
      G4Exception("G4DNASmoluchowskiReactionModel::FindReaction",
                  "G4DNASmoluchowskiReactionModel", FatalErrorInArgument,
                  exceptionDescription);
    }

    G4double __preStepSeparation =
        (__preStepPositionA - __preStepPositionB).mag();

    //===================================
    // Brownian bridge


//    if(G4Scheduler::Instance()->GetTimeStep() != __trackB.GetStep()->GetDeltaTime())
//    {
//      G4cout << G4Scheduler::Instance()->GetTimeStep() << G4endl;
//      G4cout << __trackB.GetStep()->GetDeltaTime() << G4endl;
//      assert(G4Scheduler::Instance()->GetTimeStep() == __trackB.GetStep()->GetDeltaTime());
//    }

    G4double __probabiltyOfEncounter = std::exp(
        -(__preStepSeparation - __R) * (postStepSeparation - __R) / (__D
            * (__trackB.GetStep()->GetDeltaTime())));
    G4double __selectedPOE = G4UniformRand();

    if (__selectedPOE <= __probabiltyOfEncounter) return true;
    //===================================
  }

  return false;
}
