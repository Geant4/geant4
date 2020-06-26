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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAMolecularReaction.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4VDNAReactionModel.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4MoleculeFinder.hh"
#include "G4ITReactionChange.hh"
#include "G4ITReaction.hh"

#include "G4ITTrackHolder.hh"

G4DNAMolecularReaction::G4DNAMolecularReaction()
    : G4VITReactionProcess()
    , fMolReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
    , fpReactionModel(nullptr)
{
}

G4DNAMolecularReaction::G4DNAMolecularReaction(G4VDNAReactionModel* pReactionModel)
    : G4DNAMolecularReaction()
{
    fpReactionModel = pReactionModel;
}

G4bool G4DNAMolecularReaction::TestReactibility(const G4Track &trackA,
                                                const G4Track &trackB,
                                                double currentStepTime,
                                                bool userStepTimeLimit) /*const*/
{
    const auto pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
    const auto pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();

    const G4double reactionRadius = fpReactionModel->GetReactionRadius(pMoleculeA, pMoleculeB);

    G4double separationDistance = -1.;

    if (currentStepTime == 0.)
    {
        userStepTimeLimit = false;
    }

    G4bool output = fpReactionModel->FindReaction(trackA, trackB, reactionRadius,
                                                  separationDistance, userStepTimeLimit);
    return output;
}

std::unique_ptr<G4ITReactionChange> G4DNAMolecularReaction::MakeReaction(const G4Track &trackA,
                                                                          const G4Track &trackB)
{
    std::unique_ptr<G4ITReactionChange> pChanges(new G4ITReactionChange());
    pChanges->Initialize(trackA, trackB);

    const auto pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
    const auto pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();

    const auto pReactionData = fMolReactionTable->GetReactionData(pMoleculeA, pMoleculeB);

    const G4int nbProducts = pReactionData->GetNbProducts();

    if (nbProducts)
    {
        const G4double D1 = pMoleculeA->GetDiffusionCoefficient();
        const G4double D2 = pMoleculeB->GetDiffusionCoefficient();
        const G4double sqrD1 = D1 == 0. ? 0. : std::sqrt(D1);
        const G4double sqrD2 = D2 == 0. ? 0. : std::sqrt(D2);
        const G4double inv_numerator = 1./(sqrD1 + sqrD2);
        const G4ThreeVector reactionSite = sqrD2 * inv_numerator * trackA.GetPosition()
                                         + sqrD1 * inv_numerator * trackB.GetPosition();

        for (G4int j = 0; j < nbProducts; ++j)
        {
            auto pProduct = new G4Molecule(pReactionData->GetProduct(j));
            auto pProductTrack = pProduct->BuildTrack(trackA.GetGlobalTime(), reactionSite);

            pProductTrack->SetTrackStatus(fAlive);

            G4ITTrackHolder::Instance()->Push(pProductTrack);

            pChanges->AddSecondary(pProductTrack);
            G4MoleculeFinder::Instance()->Push(pProductTrack);
        }
    }

    pChanges->KillParents(true);
    return pChanges;
}

void G4DNAMolecularReaction::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fpReactionModel = pReactionModel;
}

std::vector<std::unique_ptr<G4ITReactionChange>> G4DNAMolecularReaction::FindReaction(
    G4ITReactionSet* pReactionSet,
    const double currentStepTime,
    const double /*fGlobalTime*/,
    const bool reachedUserStepTimeLimit)
{
    std::vector<std::unique_ptr<G4ITReactionChange>> fReactionInfo;
    fReactionInfo.clear();

    if (pReactionSet == nullptr)
    {
        return fReactionInfo;
    }

    G4ITReactionPerTrackMap& reactionPerTrackMap = pReactionSet->GetReactionMap();
    for (auto tracks_i = reactionPerTrackMap.begin();
         tracks_i != reactionPerTrackMap.end();
         tracks_i = reactionPerTrackMap.begin())
    {
        G4Track* pTrackA = tracks_i->first;
        if (pTrackA->GetTrackStatus() == fStopAndKill)
        {
            continue;
        }

        G4ITReactionPerTrackPtr reactionPerTrack = tracks_i->second;
        G4ITReactionList& reactionList = reactionPerTrack->GetReactionList();

        assert(reactionList.begin() != reactionList.end());

        for (auto it = reactionList.begin(); it != reactionList.end(); it = reactionList.begin())
        {
            G4ITReactionPtr reaction(*it);
            G4Track* pTrackB = reaction->GetReactant(pTrackA);
            if (pTrackB->GetTrackStatus() == fStopAndKill)
            {
                continue;
            }

            if (pTrackB == pTrackA)
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription
                    << "The IT reaction process sent back a reaction between trackA and trackB. ";
                exceptionDescription << "The problem is trackA == trackB";
                G4Exception("G4ITModelProcessor::FindReaction",
                            "ITModelProcessor005",
                            FatalErrorInArgument,
                            exceptionDescription);
            }

            pReactionSet->SelectThisReaction(reaction);

            if (TestReactibility(*pTrackA, *pTrackB, currentStepTime, reachedUserStepTimeLimit))
            {
                auto pReactionChange = MakeReaction(*pTrackA, *pTrackB);

                if (pReactionChange)
                {
                    fReactionInfo.push_back(std::move(pReactionChange));
                    break;
                }
            }
        }
    }

    pReactionSet->CleanAllReaction();
    return fReactionInfo;
}
