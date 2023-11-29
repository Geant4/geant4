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
//
#include "G4DNASmoluchowskiReactionModel.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4UnitsTable.hh"
#include "G4Molecule.hh"
#include "G4Exp.hh"

G4DNASmoluchowskiReactionModel::G4DNASmoluchowskiReactionModel()
    : G4VDNAReactionModel()
    , fpReactionData(nullptr)
{
}

G4DNASmoluchowskiReactionModel::~G4DNASmoluchowskiReactionModel() = default;

void G4DNASmoluchowskiReactionModel::Initialise(const G4MolecularConfiguration* pMolecule,
                                                const G4Track&)
{
    fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

void G4DNASmoluchowskiReactionModel::InitialiseToPrint(const G4MolecularConfiguration* pMolecule)
{
    fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

G4double G4DNASmoluchowskiReactionModel::GetReactionRadius(const G4MolecularConfiguration* pMol1,
                                                           const G4MolecularConfiguration* pMol2)
{
    G4double __output = fpReactionTable->GetReactionData(pMol1, pMol2)->GetEffectiveReactionRadius();
    return __output;
}

G4double G4DNASmoluchowskiReactionModel::GetReactionRadius(const G4int& __i)
{
    G4double __output = (*fpReactionData)[__i]->GetEffectiveReactionRadius();
    return __output;
}

G4bool G4DNASmoluchowskiReactionModel::FindReaction(const G4Track& __trackA,
                                                    const G4Track& __trackB,
                                                    const G4double __reactionRadius,
                                                    G4double& __separationDistance,
                                                    const G4bool __alongStepReaction)
{
    const G4double R2 = __reactionRadius * __reactionRadius;
    G4double postStepSeparation = 0.;
    bool do_break = false;
    int k = 0;

    for (; k < 3; ++k)
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
        // => r^2 < R^2
        __separationDistance = std::sqrt(postStepSeparation);
        return true;
    }
    else if (__alongStepReaction == true)
    {
        //Along step check and the loop has break

        // Continue loop
        for (; k < 3; ++k)
        {
            postStepSeparation += std::pow(
                __trackA.GetPosition()[k] - __trackB.GetPosition()[k], 2);
        }
        // Use Green approach : the Brownian bridge
        __separationDistance = (postStepSeparation = std::sqrt(postStepSeparation));

        auto pMoleculeA = GetMolecule(__trackA);
        auto pMoleculeB = GetMolecule(__trackB);

        G4double D = pMoleculeA->GetDiffusionCoefficient()
                     + pMoleculeB->GetDiffusionCoefficient();

        const auto& preStepPositionA = __trackA.GetStep()->GetPreStepPoint()->GetPosition();
        const auto& preStepPositionB = __trackB.GetStep()->GetPreStepPoint()->GetPosition();

        G4double preStepSeparation = (preStepPositionA - preStepPositionB).mag();

        //===================================
        // Brownian bridge
        G4double __probabiltyOfEncounter = G4Exp(-(preStepSeparation - __reactionRadius)
                                                    * (postStepSeparation - __reactionRadius)
                                                 / (D * (__trackB.GetStep()->GetDeltaTime()))
                                            );
        G4double __selectedPOE = G4UniformRand();

        if (__selectedPOE <= __probabiltyOfEncounter)
        {
            return true;
        }
        //===================================
    }

    return false;
}
