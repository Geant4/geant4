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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAMolecularDissociation.hh"
#include "G4VITRestDiscreteProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4ParticleChange.hh"
#include "G4ITTransportationManager.hh"
#include "G4ITNavigator.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4VUserBrownianAction.hh"
//______________________________________________________________________________

G4DNAMolecularDissociation::
G4DNAMolecularDissociation(const G4String& processName,
                           G4ProcessType type)
   : G4VITRestDiscreteProcess(processName, type)
{
    // set Process Sub Type
    SetProcessSubType(fLowEnergyMolecularDecay); // DNA sub-type
    enableAlongStepDoIt = false;
    enablePostStepDoIt = true;
    enableAtRestDoIt = true;

    fVerbose = 0;

#ifdef G4VERBOSE
    if (verboseLevel > 1)
    {
        G4cout << "G4MolecularDissociationProcess constructor " << "  Name:"
               << processName << G4endl;
    }
#endif

    pParticleChange = &aParticleChange;

    fDecayAtFixedTime = true;
    fProposesTimeStep = true;
}

//______________________________________________________________________________

G4DNAMolecularDissociation::~G4DNAMolecularDissociation()
{
  delete fpBrownianAction;
}

//______________________________________________________________________________

G4bool G4DNAMolecularDissociation::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
    if (aParticleType.GetParticleType() == "Molecule")
    {
#ifdef G4VERBOSE

        if (fVerbose > 1)
        {
            G4cout << "G4MolecularDissociation::IsApplicable(";
            G4cout << aParticleType.GetParticleName() << ",";
            G4cout << aParticleType.GetParticleType() << ")" << G4endl;
        }
#endif
        return (true);
    }
    else
    {
        return false;
    }
}

//______________________________________________________________________________

G4double G4DNAMolecularDissociation::GetMeanLifeTime(const G4Track& track,
                                                     G4ForceCondition*)
{
    G4double output = GetMolecule(track)->GetDecayTime() - track.GetProperTime();
    return output > 0. ? output : 0.;
}

//______________________________________________________________________________

G4VParticleChange* G4DNAMolecularDissociation::DecayIt(const G4Track& track,
                                                       const G4Step&)
{
    aParticleChange.Initialize(track);
    auto pMotherMolecule = GetMolecule(track);
    auto pMotherMoleculeDefinition = pMotherMolecule->GetDefinition();

    if (pMotherMoleculeDefinition->GetDecayTable())
    {
        const auto pDissociationChannels = pMotherMolecule->GetDissociationChannels();

        if (pDissociationChannels == nullptr)
        {
            G4ExceptionDescription exceptionDescription;
            pMotherMolecule->PrintState();
            exceptionDescription << "No decay channel was found for the molecule : "
                                 << pMotherMolecule->GetName() << G4endl;
            G4Exception("G4DNAMolecularDissociation::DecayIt",
                        "G4DNAMolecularDissociation::NoDecayChannel",
                        FatalException,
                        exceptionDescription);
            return &aParticleChange;
        }

        auto decayVectorSize = pDissociationChannels->size();
        G4double RdmValue = G4UniformRand();

        const G4MolecularDissociationChannel* pDecayChannel = nullptr;
        size_t i = 0;
        do
        {
            pDecayChannel = (*pDissociationChannels)[i];
            if (RdmValue < pDecayChannel->GetProbability())
            {
                break;
            }
            RdmValue -= pDecayChannel->GetProbability();
            i++;
        } while (i < decayVectorSize);

        G4double decayEnergy = pDecayChannel->GetEnergy();
        auto nbProducts = pDecayChannel->GetNbProducts();

        if (decayEnergy > 0.)
        {
            aParticleChange.ProposeLocalEnergyDeposit(pDecayChannel->GetEnergy());
        }

        if (nbProducts)
        {
            std::vector<G4ThreeVector> productsDisplacement(nbProducts);
            G4ThreeVector motherMoleculeDisplacement;

            auto it = fDisplacementMap.find(pMotherMoleculeDefinition);

            if (it != fDisplacementMap.end())
            {
                auto pDisplacer = it->second.get();
                productsDisplacement = pDisplacer->GetProductsDisplacement(pDecayChannel);
                motherMoleculeDisplacement =
                        pDisplacer->GetMotherMoleculeDisplacement(pDecayChannel);
            }
            else
            {
                G4ExceptionDescription errMsg;
                errMsg << "No G4MolecularDecayProcess::theDecayDisplacementMap["
                       << pMotherMolecule->GetName() + "]";
                G4Exception("G4MolecularDecayProcess::DecayIt",
                            "DNAMolecularDecay001",
                            FatalErrorInArgument,
                            errMsg);
            }

            aParticleChange.SetNumberOfSecondaries(nbProducts);

#ifdef G4VERBOSE
            if (fVerbose)
            {
                G4cout << "Decay Process : " << pMotherMolecule->GetName()
                       << " (trackID :" << track.GetTrackID() << ") "
                       << pDecayChannel->GetName() << G4endl;
            }
#endif

            auto pNavigator = G4ITTransportationManager::GetTransportationManager()->GetNavigatorForTracking();

            for (G4int j = 0; j < nbProducts; j++)
            {
                auto pProduct = new G4Molecule(pDecayChannel->GetProduct(j));

                G4ThreeVector displacement = motherMoleculeDisplacement + productsDisplacement[j];
                double mag_displacement = displacement.mag();
                G4ThreeVector displacement_direction = displacement / (mag_displacement + 1e-30);

                double prNewSafety = DBL_MAX;

                //double step =
                pNavigator->CheckNextStep(track.GetPosition(),
                                         displacement_direction,
                                         mag_displacement,
                                         prNewSafety);

                //if(prNewSafety < mag_displacement || step < mag_displacement)
                mag_displacement = std::min(prNewSafety * 0.8, mag_displacement);

                G4ThreeVector product_pos = track.GetPosition()
                                            + displacement_direction * mag_displacement;

                //Hoang: force changing position track::
                if(fpBrownianAction != nullptr)
                {
                  fpBrownianAction->Transport(product_pos);
                }
                //Hoang: force changing position track

                const G4AffineTransform& transform = pNavigator->GetGlobalToLocalTransform();

                G4ThreeVector localPoint = transform.TransformPoint(product_pos); //track.GetPosition());
                // warning if the decayed product is outside of the volume and
                // the mother volume has no water material  (the decayed product
                // is outside of the world volume will be killed in the next step)
                if (track.GetTouchable()->GetSolid()->Inside(localPoint) !=
                    EInside::kInside)
                {
                    auto WaterMaterial = G4Material::GetMaterial("G4_WATER");
                    auto Motherlogic = track.GetTouchable()->GetVolume()->
                                       GetMotherLogical();
                    if (Motherlogic != nullptr
                       && Motherlogic->GetMaterial() != WaterMaterial)
                    {
                        G4ExceptionDescription ED;
                        ED << "The decayed product is outside of the volume : "
                           << track.GetTouchable()->GetVolume()->GetName()
                           << " with material : "<< Motherlogic->GetMaterial()
                                                       ->GetName()<< G4endl;
                        G4Exception("G4DNAMolecularDissociation::DecayIt()",
                                    "OUTSIDE_OF_MOTHER_VOLUME",
                                    JustWarning, ED);
                    }
                }

                auto pSecondary = pProduct->BuildTrack(track.GetGlobalTime(), product_pos);

                pSecondary->SetTrackStatus(fAlive);
#ifdef G4VERBOSE
                if (fVerbose)
                {
                    G4cout << "Product : " << pProduct->GetName() << G4endl;
                }
#endif
                // add the secondary track in the List
                aParticleChange.G4VParticleChange::AddSecondary(pSecondary);
            }
#ifdef G4VERBOSE
            if (fVerbose)
            {
                G4cout << "-------------" << G4endl;
            }
#endif
        }
            //DEBUG
        else if (fVerbose && decayEnergy)
        {
            G4cout << "No products for this channel" << G4endl;
            G4cout << "-------------" << G4endl;
        }
        /*
         else if(!decayEnergy && !nbProducts)
         {
         G4ExceptionDescription errMsg;
         errMsg << "There is no products and no energy specified in the molecular "
         "dissociation channel";
         G4Exception("G4MolecularDissociationProcess::DecayIt",
         "DNAMolecularDissociation002",
         FatalErrorInArgument,
         errMsg);
         }
         */
    }

    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
}

//______________________________________________________________________________

void G4DNAMolecularDissociation::SetDisplacer(Species* pSpecies, Displacer* pDisplacer)
{
    fDisplacementMap.emplace(pSpecies, std::unique_ptr<Displacer>(pDisplacer));
}

//______________________________________________________________________________

G4DNAMolecularDissociation::Displacer* G4DNAMolecularDissociation::GetDisplacer(Species* pSpecies)
{
    return fDisplacementMap[pSpecies].get();
}

//______________________________________________________________________________

G4double G4DNAMolecularDissociation::PostStepGetPhysicalInteractionLength(const G4Track&,
                                                                          G4double,
                                                                          G4ForceCondition*)
{
    return 0; //-1*(-log(1-G4UniformRand())*100*1e-15*s);
}

//______________________________________________________________________________

void G4DNAMolecularDissociation::SetVerbose(G4int verbose)
{
    fVerbose = verbose;
}

//______________________________________________________________________________

G4VParticleChange* G4DNAMolecularDissociation::AtRestDoIt(const G4Track& track,
                                                          const G4Step& step)
{
    ClearNumberOfInteractionLengthLeft();
    ClearInteractionTimeLeft();
    return DecayIt(track, step);
}

//______________________________________________________________________________

G4double G4DNAMolecularDissociation::AtRestGetPhysicalInteractionLength(const G4Track& track,
                                                                        G4ForceCondition* condition)
{
    if (fDecayAtFixedTime)
    {
        return GetMeanLifeTime(track, condition);
    }

    return G4VITRestDiscreteProcess::AtRestGetPhysicalInteractionLength(track, condition);
}

//______________________________________________________________________________

G4VParticleChange* G4DNAMolecularDissociation::PostStepDoIt(const G4Track& track,
                                                            const G4Step& step)
{
    return AtRestDoIt(track, step);
}

//______________________________________________________________________________

G4double G4DNAMolecularDissociation::GetMeanFreePath(const G4Track&,
                                                     G4double,
                                                     G4ForceCondition*)
{
    return 0;
}
