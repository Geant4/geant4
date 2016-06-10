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
#include "G4MoleculeCounter.hh"

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
  return 1 + std::exp(std::log(10.)*
                 (Y(density) *
                     (C(temperature) + (S(temperature) - 1)*std::log(density)/std::log(10.))
                     + D(temperature) + std::log(density)/std::log(10.)));
}

//------------------------------------------------------------------------------

G4DNAElectronHoleRecombination::G4DNAElectronHoleRecombination() :
    G4VITRestDiscreteProcess("G4DNAElectronHoleRecombination", fElectromagnetic)
{
  Create();
//  G4cout << epsilon(1.0095, 298.) << G4endl;
//  G4cout << epsilon(1., 293.15) << G4endl;
//  G4cout << epsilon(0.9277, 423.) << G4endl;
//  G4cout << epsilon(0.816, 523.) << G4endl;
//  G4cout << epsilon(0.6, 623) << G4endl;
}

G4DNAElectronHoleRecombination::~G4DNAElectronHoleRecombination()
{
}

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
  fpMoleculeDensity = 0;

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

void G4DNAElectronHoleRecombination::StartTracking(G4Track* track)
{
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState.reset(new State());
  G4VITProcess::StartTracking(track);
}

//______________________________________________________________________________

void G4DNAElectronHoleRecombination::MakeReaction(const G4Track& track)
{
  fParticleChange.Initialize(track);
  State* state = fpState->GetState<State>();
  double random = state->fSampleProba;
  std::vector<ReactionProfile>& reactants = state->fReactants;

  G4Track* selected_reactant = 0;

  for(size_t i = 0; i < reactants.size(); ++i)
  {
    if(reactants[i].fElectron->GetTrackStatus() != fAlive) continue;
    if(reactants[i].fProbability > random)
    {
      selected_reactant = reactants[i].fElectron;
    }
    break;
  }

  //  G4cout << "MakeReaction with charge ="
  //       << GetMolecule(track)->GetCharge() << G4endl;

  if(selected_reactant)
  {
    //    G4cout << " Will react with TID = " << selected_reactant->GetTrackID()
    //           << G4endl;

    if(G4MoleculeCounter::InUse())
      G4MoleculeCounter::Instance()->
        RemoveAMoleculeAtTime(GetMolecule(track)->GetMolecularConfiguration(),
                              track.GetGlobalTime());
    GetMolecule(track)->ChangeConfigurationToLabel("H2Ovib");

    if(G4MoleculeCounter::InUse())
      G4MoleculeCounter::Instance()->
        AddAMoleculeAtTime(GetMolecule(track)->GetMolecularConfiguration(),
                           track.GetGlobalTime());

    //  fParticleChange.ProposeTrackStatus(fStopAndKill);
    fParticleChange.ProposeTrackStatus(fStopButAlive);

    selected_reactant->SetTrackStatus(fStopAndKill);
    //  G4TrackList::Pop(selected_reactant);
    //  G4ITTrackHolder::Instance()->PushToKill(selected_reactant);

  }
  else
  {
    fParticleChange.ProposeTrackStatus(fStopButAlive);
  }
}

//______________________________________________________________________________

G4bool G4DNAElectronHoleRecombination::FindReactant(const G4Track& track)
{
  if(GetMolecule(track)->GetCharge() <= 0)
  {
    // G4cout << "La charge est negative ou nulle !! " << G4endl;
    return false;
  }

  const std::vector<double>* densityTable =
      G4DNAMolecularMaterial::Instance()->GetDensityTableFor(track.GetMaterial());

  double temperature = track.GetMaterial()->GetTemperature();
  double density = (*densityTable)[track.GetMaterial()->GetIndex()] /
                    ( g/(1e-2*m*1e-2*m*1e-2*m) );
  double eps = epsilon(density, temperature);

  // G4cout << " temperature = " << temperature << G4endl;
  // G4cout << " density = " << density << G4endl;
  // G4cout << " eps = " << eps << G4endl;

  double onsager_radius = onsager_constant * 1. / (temperature * eps);

  G4Molecule e_aq(G4Electron_aq::Definition());

  G4KDTreeResultHandle results = G4MoleculeFinder::Instance()
      ->FindNearestInRange(track.GetPosition(),
                           e_aq.GetMoleculeID(),
                           10. * onsager_radius);

  // double distance = -1.;
  // double probability = -1.;

  if(results == 0 || results->GetSize() == 0)
  {
    //    G4cout << "rien trouve a moins de 10 rc" << G4endl;
    return false;
  }

  results->Sort();

  State* state = fpState->GetState<State>();
  std::vector<ReactionProfile>& reactants = state->fReactants;
  state->fSampleProba = G4UniformRand();

  reactants.resize(results->GetSize());

  for(size_t i = 0; results->End() == false; results->Next(), ++i)
  {
    reactants[i].fElectron = results->GetItem<G4IT>()->GetTrack();
    reactants[i].fDistance = std::sqrt(results->GetDistanceSqr());

    if(reactants[i].fDistance != 0)
    {
      reactants[i].fProbability = 1.
          - std::exp(-onsager_radius / reactants[i].fDistance);
    }
    else
    {
      reactants[i].fProbability = 1.;
    }

    //  G4cout << "dis = "
    //         << reactants[i].fDistance << " prob = " << reactants[i].fProbability << G4endl;
  }

  if(results->GetSize() != 0 && reactants.empty())
  {
    G4cout << "Size is = " << results->GetSize() << G4endl;
    abort();
  }

  if(reactants.empty()) return false;

  //  G4cout << " reactants[0].fDistance =" << reactants[0].fDistance
  //        << " onsager_radius = " << onsager_radius << "\t";
  //
  //  G4cout << " reactants[0].fProbability =" << reactants[0].fProbability
  //        << "state->fSampleProba = " << state->fSampleProba << "\t";

  if(reactants[0].fProbability > state->fSampleProba) return true;
  return false;
}

//______________________________________________________________________________

G4bool
G4DNAElectronHoleRecombination::
IsApplicable(const G4ParticleDefinition& particle)
{
  if(G4Threading::IsMasterThread())
  {
    G4MoleculeDefinition* H2O = G4MoleculeTable::Instance()
        ->GetMoleculeDefinition("H2O", false);

    if(H2O) // if this condition does not hold => process cannot be applied
    {
      G4MolecularConfiguration* vib =
          G4H2O::Definition()->NewConfiguration("H2Ovib");

      assert(vib != 0);

      G4MolecularConfiguration* H2 = G4MoleculeTable::Instance()
          ->GetConfiguration("H2", false);
      G4MolecularConfiguration* OH = G4MoleculeTable::Instance()
          ->GetConfiguration("OH", false);
      G4MolecularConfiguration* H = G4MoleculeTable::Instance()
          ->GetConfiguration("H", false);

      double probaRemaining = 1.;

      if(OH || H2)
      {
        G4MolecularDissociationChannel* diss2 =
                  new G4MolecularDissociationChannel("H2Ovib_DissociativeDecay1");
        if(H2)
        {
          diss2->AddProduct(H2);
        }
        if(OH)
        {
          diss2->AddProduct(OH);
          diss2->AddProduct(OH);
        }

        double proba = 0.15;
        diss2->SetProbability(proba);
        probaRemaining -= proba;
        diss2->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                                   B1A1_DissociationDecay);
        G4H2O::Definition()->AddDecayChannel(vib, diss2);
      }

      if(OH || H)
      {
        G4MolecularDissociationChannel* diss3 =
            new G4MolecularDissociationChannel("H2Ovib_DissociativeDecay2");
        if(OH)
        {
          diss3->AddProduct(OH);
        }
        if(H)
        {
          diss3->AddProduct(H);
        }
        double proba = 0.55;
        diss3->SetProbability(proba);
        probaRemaining -= proba;
        diss3->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                                   A1B1_DissociationDecay);
        G4H2O::Definition()->AddDecayChannel(vib, diss3);
      }

      G4MolecularDissociationChannel* channel1 =
          new G4MolecularDissociationChannel("H2Ovib_NonDissociative");
      channel1->SetProbability(probaRemaining);
      G4H2O::Definition()->AddDecayChannel(vib, channel1);
    }
  }

  if(particle.GetParticleName() == "H2O") return true;
  return false;
}

//______________________________________________________________________________

G4double G4DNAElectronHoleRecombination::GetMeanFreePath(const G4Track& track,
                                                         G4double,
                                                         G4ForceCondition*)
{
  //  G4cout << "track ID = " << track.GetTrackID()
  //      <<  " G4DNAElectronHoleRecombination --> ";
    if(FindReactant(track))
    {
  //    G4cout << " SELECTED " ;
  //    G4cout << G4endl;
      return 0;
    }

  //  G4cout << " NOT SELECTED " << G4endl;
    return DBL_MAX;
}

//______________________________________________________________________________

G4double G4DNAElectronHoleRecombination::GetMeanLifeTime(const G4Track& track,
                                                         G4ForceCondition*)
{
//  G4cout << "track ID = " << track.GetTrackID()
//      <<  " G4DNAElectronHoleRecombination --> ";
  if(FindReactant(track))
  {
//    G4cout << " SELECTED " ;
//    G4cout << G4endl;
    return 0;
  }

//  G4cout << " NOT SELECTED " << G4endl;
  return DBL_MAX;
}

G4VParticleChange* G4DNAElectronHoleRecombination::PostStepDoIt(const G4Track& track,
                                                                const G4Step& step)
{
  return AtRestDoIt(track, step);
}
