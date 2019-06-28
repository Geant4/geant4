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
/*
 * G4DNAElectronHoleRecombination.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: mkaramit
 */

#include <G4DNAElectronHoleRecombination.hh>
#include <G4MoleculeFinder.hh>
#include "G4PhysicalConstants.hh"
#include "G4Electron_aq.hh"
#include "G4MoleculeTable.hh"
#include "G4MolecularDissociationChannel.hh"
#include "G4H2.hh"
#include "G4H2O.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4SystemOfUnits.hh"
#include "G4VMoleculeCounter.hh"
#include "G4Exp.hh"

static double onsager_constant = e_squared / (4. * pi * epsilon0 * k_Boltzmann);

//------------------------------------------------------------------------------

// Parameterisation of dielectric constant vs temperature and density

double Y(double density)
{
    return 1. / (1. + 0.0012 / (density * density));
}

double A(double temperature)
{
    double temp_inverse = 1 / temperature;
    return 0.7017
           + 642.0 * temp_inverse
           - 1.167e5 * temp_inverse * temp_inverse
           + 9.190e6 * temp_inverse * temp_inverse * temp_inverse;
}

double B(double temperature)
{
    double temp_inverse = 1 / temperature;
    return -2.71
           + 275.4 * temp_inverse
           + 0.3245e5 * temp_inverse * temp_inverse;
}

double S(double temp)
{
    double temp_inverse = 1 / temp;

    return 1.667
           - 11.41 * temp_inverse
           - 35260.0 * temp_inverse * temp_inverse;
}

double C(double temp)
{
    return A(temp) - B(temp) - 3;
}

double D(double temp)
{
    return B(temp) + 3;
}

double epsilon(double density, double temperature)
{
    return 1 + G4Exp(std::log(10.) *
                     (Y(density) *
                      (C(temperature) + (S(temperature) - 1) * std::log(density) / std::log(10.))
                      + D(temperature) + std::log(density) / std::log(10.)));
}

//------------------------------------------------------------------------------

G4DNAElectronHoleRecombination::G4DNAElectronHoleRecombination()
        : G4VITRestDiscreteProcess("G4DNAElectronHoleRecombination",
                                   fElectromagnetic)
{
    Create();
}

G4DNAElectronHoleRecombination::~G4DNAElectronHoleRecombination() = default;

void G4DNAElectronHoleRecombination::Create()
{
    pParticleChange = &fParticleChange;
    enableAtRestDoIt = true;
    enableAlongStepDoIt = false;
    enablePostStepDoIt = true;

    SetProcessSubType(60);

    G4VITProcess::SetInstantiateProcessState(false);
    // ie G4DNAElectronHoleRecombination uses a state class
    // inheriting from G4ProcessState

    fIsInitialized = false;
    fProposesTimeStep = true;
    fpMoleculeDensity = nullptr;

    verboseLevel = 0;
}

//______________________________________________________________________________

G4VParticleChange*
G4DNAElectronHoleRecombination::AtRestDoIt(const G4Track& track,
                                           const G4Step& /*stepData*/)
{
    fParticleChange.Initialize(track);
    ClearInteractionTimeLeft();
    ClearNumberOfInteractionLengthLeft();
    MakeReaction(track);
    return &fParticleChange;
}

//______________________________________________________________________________

void G4DNAElectronHoleRecombination::StartTracking(G4Track* pTrack)
{
    G4VProcess::StartTracking(pTrack);
    G4VITProcess::fpState.reset(new State());
    G4VITProcess::StartTracking(pTrack);
}

//______________________________________________________________________________

void G4DNAElectronHoleRecombination::MakeReaction(const G4Track& track)
{
    fParticleChange.Initialize(track);
    auto pState = fpState->GetState<State>();
    double random = pState->fSampleProba;
    std::vector<ReactantInfo>& reactants = pState->fReactants;

    G4Track* pSelectedReactant = nullptr;

    for (const auto& reactantInfo : reactants)
    {
        if (reactantInfo.fElectron->GetTrackStatus() != fAlive)
        {
            continue;
        }
        if (reactantInfo.fProbability > random)
        {
            pSelectedReactant = reactantInfo.fElectron;
        }
        break;
    }

    if (pSelectedReactant)
    {
        if (G4VMoleculeCounter::InUse())
        {
            G4VMoleculeCounter::Instance()->
                    RemoveAMoleculeAtTime(GetMolecule(track)->GetMolecularConfiguration(),
                                          track.GetGlobalTime(),
                                          &(track.GetPosition()));
        }
        GetMolecule(track)->ChangeConfigurationToLabel("H2Ovib");

        if (G4VMoleculeCounter::InUse())
        {
            G4VMoleculeCounter::Instance()->
                    AddAMoleculeAtTime(GetMolecule(track)->GetMolecularConfiguration(),
                                       track.GetGlobalTime(),
                                       &(track.GetPosition()));
        }

        //  fParticleChange.ProposeTrackStatus(fStopAndKill);
        fParticleChange.ProposeTrackStatus(fStopButAlive);

        pSelectedReactant->SetTrackStatus(fStopAndKill);
        //  G4TrackList::Pop(pSelectedReactant);
        //  G4ITTrackHolder::Instance()->PushToKill(pSelectedReactant);
    }
    else
    {
        fParticleChange.ProposeTrackStatus(fStopButAlive);
    }
}

//______________________________________________________________________________

G4bool G4DNAElectronHoleRecombination::FindReactant(const G4Track& track)
{
    if (GetMolecule(track)->GetCharge() <= 0)
    {
        // cannot react with electrons
        return false;
    }

    const auto pDensityTable =
            G4DNAMolecularMaterial::Instance()->GetDensityTableFor(track.GetMaterial());

    double temperature = track.GetMaterial()->GetTemperature();
    double density = (*pDensityTable)[track.GetMaterial()->GetIndex()] / (g / (1e-2 * m * 1e-2 * m * 1e-2 * m));
    double eps = epsilon(density, temperature);

    double onsagerRadius = onsager_constant * 1. / (temperature * eps);

    G4Molecule e_aq(G4Electron_aq::Definition());

    G4KDTreeResultHandle results = G4MoleculeFinder::Instance()
            ->FindNearestInRange(track.GetPosition(),
                                 e_aq.GetMoleculeID(),
                                 10. * onsagerRadius);

    if (results == 0 || results->GetSize() == 0)
    {
        return false;
    }

    results->Sort();

    auto pState = fpState->GetState<State>();
    std::vector<ReactantInfo>& reactants = pState->fReactants;
    pState->fSampleProba = G4UniformRand();

    reactants.resize(results->GetSize());

    for (size_t i = 0; !results->End(); results->Next(), ++i)
    {
        reactants[i].fElectron = results->GetItem<G4IT>()->GetTrack();
        reactants[i].fDistance = std::sqrt(results->GetDistanceSqr());

        if (reactants[i].fDistance != 0)
        {
            reactants[i].fProbability = 1. - G4Exp(-onsagerRadius / reactants[i].fDistance);
        }
        else
        {
            reactants[i].fProbability = 1.;
        }
    }

    return reactants.empty() ? false : reactants[0].fProbability > pState->fSampleProba;
}

//______________________________________________________________________________

void G4DNAElectronHoleRecombination::BuildDissociationChannels()
{
    auto pMoleculeTable = G4MoleculeTable::Instance();

    auto pH2ODef = pMoleculeTable->GetMoleculeDefinition("H2O", false);

    if (pH2ODef == nullptr) // if this condition does not hold => process cannot be applied
    {
        return;
    }

    const auto pH2Ovib = G4H2O::Definition()->NewConfiguration("H2Ovib");

    assert(pH2Ovib != nullptr);

    const auto pH2 = pMoleculeTable->GetConfiguration("H2", false);
    const auto pOH = pMoleculeTable->GetConfiguration("OH", false);
    const auto pH = pMoleculeTable->GetConfiguration("H", false);

    double probaRemaining = 1.;

    if (pOH || pH2)
    {
        auto pDissocationChannel2 =
                new G4MolecularDissociationChannel("H2Ovib_DissociativeDecay1");
        if (pH2)
        {
            pDissocationChannel2->AddProduct(pH2);
        }
        if (pOH)
        {
            pDissocationChannel2->AddProduct(pOH);
            pDissocationChannel2->AddProduct(pOH);
        }

        double channelProbability = 0.15;
        pDissocationChannel2->SetProbability(channelProbability);
        probaRemaining -= channelProbability;
        pDissocationChannel2->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                                                  B1A1_DissociationDecay);
        pH2ODef->AddDecayChannel(pH2Ovib, pDissocationChannel2);
    }

    if (pOH || pH)
    {
        auto pDissociationChannel3 =
                new G4MolecularDissociationChannel("H2Ovib_DissociativeDecay2");
        if (pOH)
        {
            pDissociationChannel3->AddProduct(pOH);
        }
        if (pH)
        {
            pDissociationChannel3->AddProduct(pH);
        }
        double channelProbability = 0.55;
        pDissociationChannel3->SetProbability(channelProbability);
        probaRemaining -= channelProbability;
        pDissociationChannel3->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                                                   A1B1_DissociationDecay);
        pH2ODef->AddDecayChannel(pH2Ovib, pDissociationChannel3);
    }

    auto pDissociationChannel1 =
            new G4MolecularDissociationChannel("H2Ovib_NonDissociative");
    pDissociationChannel1->SetProbability(probaRemaining);
    pH2ODef->AddDecayChannel(pH2Ovib, pDissociationChannel1);
}

G4bool
G4DNAElectronHoleRecombination::
IsApplicable(const G4ParticleDefinition& particle)
{
    if (&particle != G4H2O::DefinitionIfExists())
    {
        return false;
    }

    if (G4Threading::IsMasterThread())
    {
        BuildDissociationChannels();
    }

    return true;
}

//______________________________________________________________________________

G4double G4DNAElectronHoleRecombination::GetMeanFreePath(const G4Track& track,
                                                         G4double,
                                                         G4ForceCondition*)
{
    if (FindReactant(track))
    {
        return 0.;
    }

    return DBL_MAX;
}

//______________________________________________________________________________

G4double G4DNAElectronHoleRecombination::GetMeanLifeTime(const G4Track& track,
                                                         G4ForceCondition*)
{
    if (FindReactant(track))
    {
        return 0.;
    }

    return DBL_MAX;
}

G4VParticleChange* G4DNAElectronHoleRecombination::PostStepDoIt(const G4Track& track,
                                                                const G4Step& step)
{
    return AtRestDoIt(track, step);
}
