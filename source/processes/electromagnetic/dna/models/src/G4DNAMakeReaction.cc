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


#include "G4DNAMakeReaction.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4VDNAReactionModel.hh"
#include "G4Molecule.hh"
#include "G4MoleculeFinder.hh"
#include "G4ITReactionChange.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4ITReaction.hh"
#include "G4DNAIndependentReactionTimeStepper.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"

G4DNAMakeReaction::G4DNAMakeReaction()
    : G4VITReactionProcess()
    , fMolReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
    , fpReactionModel(nullptr)
    , fpTimeStepper(nullptr)
    , fTimeStep(0)
{
}

G4DNAMakeReaction::G4DNAMakeReaction(G4VDNAReactionModel* pReactionModel)
        : G4DNAMakeReaction()
{
    fpReactionModel = pReactionModel;
}

void G4DNAMakeReaction::SetTimeStepComputer(G4VITTimeStepComputer* pStepper)
{
    fpTimeStepper = pStepper;
}

G4bool G4DNAMakeReaction::TestReactibility(const G4Track& /*trackA*/,
                                           const G4Track& /*trackB*/,
                                           G4double currentStepTime,
                                           G4bool /*userStepTimeLimit*/) /*const*/
{
    fTimeStep = currentStepTime;
    return true;
}

std::unique_ptr<G4ITReactionChange>
G4DNAMakeReaction::MakeReaction(const G4Track &trackA,
                                const G4Track &trackB)
{
    G4Track& tA = const_cast<G4Track&>(trackA);
    G4Track& tB = const_cast<G4Track&>(trackB);
    UpdatePositionForReaction( tA , tB );//TODO: should change it

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
        const G4ThreeVector reactionSite = sqrD2 * inv_numerator * tA.GetPosition()
                                         + sqrD1 * inv_numerator * tB.GetPosition();

        G4double u = G4UniformRand();
        auto randP = (1-u) * tA.GetPosition() + u * tB.GetPosition();

        for (G4int j = 0; j < nbProducts; ++j)
        {
            auto pProduct = new G4Molecule(pReactionData->GetProduct(j));
            auto pProductTrack = pProduct->BuildTrack(trackA.GetGlobalTime(), (reactionSite + randP)/2);
            pProductTrack->SetTrackStatus(fAlive);
            G4ITTrackHolder::Instance()->Push(pProductTrack);
            pChanges->AddSecondary(pProductTrack);
        }
    }
    pChanges->KillParents(true);
    return pChanges;
}

void G4DNAMakeReaction::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fpReactionModel = pReactionModel;
}

void G4DNAMakeReaction::UpdatePositionForReaction(G4Track& trackA,
                                           G4Track& trackB)
{
    const auto pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
    const auto pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();
    G4double D1 = pMoleculeA->GetDiffusionCoefficient();
    G4double D2 = pMoleculeB->GetDiffusionCoefficient();

    G4double reactionRadius = fpReactionModel->GetReactionRadius( pMoleculeA, pMoleculeB );
    G4ThreeVector p1 = trackA.GetPosition();
    G4ThreeVector p2 = trackB.GetPosition();

    G4ThreeVector S1 = p1 - p2;
    G4double distance = S1.mag();

    if(D1 == 0)
    {
        trackB.SetPosition(p1);
        return;
    }
    else if(D2 == 0)
    {
        trackA.SetPosition(p2);
        return;
    }

    if(distance == 0)
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription << "Two particles are overlap: "
                             <<GetMolecule(trackA)->GetName()
                             <<" and "<<GetMolecule(trackB)->GetName()
                             <<" at "<<trackA.GetPosition();
        G4Exception("G4DNAMakeReaction::PrepareForReaction()",
                    "G4DNAMakeReaction003",
                    FatalErrorInArgument,exceptionDescription);
    }
    S1.setMag(reactionRadius);

    const G4double dt = fTimeStep;//irt - actualize molecule time

    if(dt > 0)// irt > 0
    {
        G4double s12 = 2.0 * D1 * dt;
        G4double s22 = 2.0 * D2 * dt;
        G4double sigma = s12 + ( s12 * s12 ) / s22;
        G4double alpha = reactionRadius * distance / (2 * (D1 + D2) * dt );

        G4ThreeVector S2 = (p1 + ( s12 / s22 ) * p2) +
                           G4ThreeVector(G4RandGauss::shoot(0.0, sigma),
                                         G4RandGauss::shoot(0.0, sigma),
                                         G4RandGauss::shoot(0.0, sigma));

        S1.setPhi(rad * G4UniformRand() * 2.0 * CLHEP::pi);

        S1.setTheta(rad * std::acos( 1.0 + (1. / alpha) *
                                           std::log(1.0 - G4UniformRand() *
                                                          (1.-std::exp(-2.0 * alpha)))));

        const G4ThreeVector R1 = (D1 * S1 + D2 * S2) / (D1 + D2);
        const G4ThreeVector R2 = D2 * (S2 - S1) / (D1 + D2);

        trackA.SetPosition(R1);
        trackB.SetPosition(R2);
    }
}

std::vector<std::unique_ptr<G4ITReactionChange>>
G4DNAMakeReaction::FindReaction(G4ITReactionSet* pReactionSet,
                                const G4double currentStepTime,
                                const G4double /*globalTime*/,
                                const G4bool /*reachedUserStepTimeLimit*/)
{
    std::vector<std::unique_ptr<G4ITReactionChange>> ReactionInfo;
    ReactionInfo.clear();
    auto stepper = dynamic_cast<G4DNAIndependentReactionTimeStepper*>(fpTimeStepper);
    if(stepper != nullptr){
      auto pReactionChange = stepper->
                             FindReaction(pReactionSet,currentStepTime);
      if (pReactionChange != nullptr)
      {
        ReactionInfo.push_back(std::move(pReactionChange));
      }
    }
    return ReactionInfo;
}
