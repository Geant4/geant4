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

#include <iomanip>

#include "G4DNAMolecularReactionTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"
#include "G4VDNAReactionModel.hh"
#include "G4MoleculeHandleManager.hh"
#include "G4MoleculeTable.hh"
#include "G4MolecularConfiguration.hh"
#include "G4ReactionTableMessenger.hh"
#include "G4IosFlagsSaver.hh"
#include "G4Exp.hh"

using namespace std;

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::fpInstance(nullptr);

G4DNAMolecularReactionData::G4DNAMolecularReactionData()
    : fpReactant1(nullptr)
    , fpReactant2(nullptr)
    , fObservedReactionRate(0.)
    , fActivationRate(0.)
    , fDiffusionRate(0.)
    , fOnsagerRadius(0.)
    , fReactionRadius(0.)
    , fEffectiveReactionRadius(0.)
    , fProbability(0.)
    , fType(0)
    , fReactionID(0)
{
}

G4DNAMolecularReactionData::G4DNAMolecularReactionData(G4double reactionRate,
                                                       Reactant* pReactant1,
                                                       Reactant* pReactant2)
    : fpReactant1(pReactant1)
    , fpReactant2(pReactant2)
    , fObservedReactionRate(reactionRate)
    , fActivationRate(0.)
    , fDiffusionRate(0.)
    , fOnsagerRadius(0.)
    , fReactionRadius(0.)
    , fEffectiveReactionRadius(0.)
    , fProbability(0.)
    , fType(0)
    , fReactionID(0)
{
    ComputeEffectiveRadius();
}

G4DNAMolecularReactionData::G4DNAMolecularReactionData(G4double reactionRate,
                                                       const G4String& reactant1,
                                                       const G4String& reactant2)
    : fpReactant1(nullptr)
    , fpReactant2(nullptr)
    , fObservedReactionRate(reactionRate)
    , fActivationRate(0.)
    , fDiffusionRate(0.)
    , fOnsagerRadius(0.)
    , fReactionRadius(0.)
    , fEffectiveReactionRadius(0.)
    , fProbability(0.)
    , fType(0)
    , fReactionID(0)
{
    SetReactant1(reactant1);
    SetReactant2(reactant2);
    ComputeEffectiveRadius();
}

G4DNAMolecularReactionData::~G4DNAMolecularReactionData()
{
     fProducts.clear();
}

void G4DNAMolecularReactionData::ComputeEffectiveRadius()
{
    G4double sumDiffCoeff = 0.;

    if (fpReactant1 == fpReactant2)
    {
        sumDiffCoeff = fpReactant1->GetDiffusionCoefficient();
        fEffectiveReactionRadius = fObservedReactionRate / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
    }
    else
    {
        sumDiffCoeff = fpReactant1->GetDiffusionCoefficient()
                     + fpReactant2->GetDiffusionCoefficient();
        fEffectiveReactionRadius = fObservedReactionRate / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
    }

    fReactionID = 0;
    fReactionRadius = fEffectiveReactionRadius;
    fOnsagerRadius = (fpReactant1->GetCharge() * fpReactant2->GetCharge())/(4*pi*epsilon0*k_Boltzmann) / (293.15 * 80.1) ;
    fProbability = 1;

}

int G4DNAMolecularReactionData::GetReactionID() const
{
    return fReactionID;
}

void G4DNAMolecularReactionData::SetReactionID(int ID)
{
    fReactionID = ID;
}

void G4DNAMolecularReactionData::SetReactant1(Reactant* pReactive)
{
    fpReactant1 = pReactive;
}

void G4DNAMolecularReactionData::SetReactant2(Reactant* pReactive)
{
    fpReactant2 = pReactive;
}

void G4DNAMolecularReactionData::SetReactants(Reactant* pReactant1,
                                              Reactant* pReactant2)
{
    fpReactant1 = pReactant1;
    fpReactant2 = pReactant2;
}

void G4DNAMolecularReactionData::AddProduct(Reactant* pMolecule)
{
    fProducts.push_back(pMolecule);
}

G4int G4DNAMolecularReactionData::GetNbProducts() const
{
    return (G4int)fProducts.size();
}

G4DNAMolecularReactionData::Reactant* G4DNAMolecularReactionData::GetProduct(G4int i) const
{
    return fProducts[i];
}

const G4DNAMolecularReactionData::ReactionProducts* G4DNAMolecularReactionData::GetProducts() const
{
    return &fProducts;
}

void G4DNAMolecularReactionData::RemoveProducts()
{
    fProducts.clear();
}

void G4DNAMolecularReactionData::SetReactant1(const G4String& reactive)
{
    fpReactant1 = G4MoleculeTable::Instance()->GetConfiguration(reactive);
}
void G4DNAMolecularReactionData::SetReactant2(const G4String& reactive)
{
    fpReactant2 = G4MoleculeTable::Instance()->GetConfiguration(reactive);
}
void G4DNAMolecularReactionData::SetReactants(const G4String& reactant1,
                                              const G4String& reactant2)
{
    fpReactant1 = G4MoleculeTable::Instance()->GetConfiguration(reactant1);
    fpReactant2 = G4MoleculeTable::Instance()->GetConfiguration(reactant2);
}

G4DNAMolecularReactionData::ReactantPair G4DNAMolecularReactionData::GetReactants()
{
    return std::make_pair(fpReactant1, fpReactant2);
}

G4DNAMolecularReactionData::Reactant* G4DNAMolecularReactionData::GetReactant1() const
{
    return fpReactant1;
}

G4DNAMolecularReactionData::Reactant* G4DNAMolecularReactionData::GetReactant2() const
{
    return fpReactant2;
}

void G4DNAMolecularReactionData::SetObservedReactionRateConstant(G4double rate)
{
    fObservedReactionRate = rate;
}

G4double G4DNAMolecularReactionData::GetObservedReactionRateConstant() const
{
    return fObservedReactionRate;
}

G4double G4DNAMolecularReactionData::GetActivationRateConstant() const
{
    return fActivationRate;
}

G4double G4DNAMolecularReactionData::GetDiffusionRateConstant() const
{
    return fDiffusionRate;
}

void G4DNAMolecularReactionData::SetReactionRadius(G4double radius)
{
    fReactionRadius = radius;
    fEffectiveReactionRadius = -fOnsagerRadius / (1-exp(fOnsagerRadius / fReactionRadius));
}

G4double G4DNAMolecularReactionData::GetReactionRadius() const
{
    return fReactionRadius;
}

void G4DNAMolecularReactionData::SetEffectiveReactionRadius(G4double radius)
{
    fEffectiveReactionRadius = radius;
}

G4double G4DNAMolecularReactionData::GetEffectiveReactionRadius() const
{
    return fEffectiveReactionRadius;
}

G4double G4DNAMolecularReactionData::GetOnsagerRadius() const
{
    return fOnsagerRadius;
}

G4double G4DNAMolecularReactionData::GetProbability() const
{
    return fProbability;
}

void G4DNAMolecularReactionData::SetProbability(G4double prob)
{
    fProbability = prob;
}

void G4DNAMolecularReactionData::SetReactionType(G4int type)
{
    G4double sumDiffCoeff = 0.;

    if(type == 1)
    {

        sumDiffCoeff = fpReactant1->GetDiffusionCoefficient() + 
                       fpReactant2->GetDiffusionCoefficient();

        fReactionRadius = fpReactant1->GetVanDerVaalsRadius() +
                          fpReactant2->GetVanDerVaalsRadius();

        G4double Rs = 0.29 * nm;
        if(fOnsagerRadius == 0) // Type II
        {
            fEffectiveReactionRadius = fReactionRadius;
            fDiffusionRate = 4 * pi * sumDiffCoeff * fReactionRadius * Avogadro;
            if (fpReactant1 == fpReactant2) fDiffusionRate/=2;
            fActivationRate = fDiffusionRate * fObservedReactionRate / (fDiffusionRate - fObservedReactionRate);
            fProbability =  Rs / (Rs + (fDiffusionRate / fActivationRate) * (fReactionRadius + Rs));

        }else{ // Type IV
            fEffectiveReactionRadius = -fOnsagerRadius/(1-exp(fOnsagerRadius/fReactionRadius));
            fDiffusionRate = 4 * pi * sumDiffCoeff * fEffectiveReactionRadius * Avogadro;
            if (fpReactant1 == fpReactant2) fDiffusionRate/=2;

            fActivationRate = fDiffusionRate * fObservedReactionRate / (fDiffusionRate - fObservedReactionRate);
            fProbability = Rs / (Rs + (fDiffusionRate / fActivationRate) * (fEffectiveReactionRadius + Rs));
        }
    }

    fType = type;
}

G4int G4DNAMolecularReactionData::GetReactionType() const
{
    return fType;
}

void G4DNAMolecularReactionData::AddProduct(const G4String& molecule)
{
    fProducts.push_back(G4MoleculeTable::Instance()->GetConfiguration(molecule));
}

double G4DNAMolecularReactionData::PolynomialParam(double temp_K, std::vector<double> P)
{
    double inv_temp = 1. / temp_K;

    return pow(10,
               P[0] + P[1] * inv_temp + P[2] * pow(inv_temp, 2)
               + P[3] * pow(inv_temp, 3) + P[4] * pow(inv_temp, 4))
        * (1e-3 * CLHEP::m3 / (CLHEP::mole * CLHEP::s));
}

double G4DNAMolecularReactionData::ArrehniusParam(double temp_K, std::vector<double> P)
{
    return P[0] * G4Exp(P[1] / temp_K)*
        (1e-3 * CLHEP::m3 / (CLHEP::mole * CLHEP::s));
}

double G4DNAMolecularReactionData::ScaledParameterization(double temp_K,
                                                          double temp_init,
                                                          double rateCste_init)
{
    double D0 = G4MolecularConfiguration::DiffCoeffWater(temp_init);
    double Df = G4MolecularConfiguration::DiffCoeffWater(temp_K);
    return Df * rateCste_init / D0;
}

//==============================================================================
// REACTION TABLE
//==============================================================================

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::GetReactionTable()
{
    if (!fpInstance)
    {
        fpInstance = new G4DNAMolecularReactionTable();
    }
    return fpInstance;
}

//_____________________________________________________________________________________

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::Instance()
{
    if (!fpInstance)
    {
        fpInstance = new G4DNAMolecularReactionTable();
    }
    return fpInstance;
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::DeleteInstance()
{
    if (fpInstance)
    {
        delete fpInstance;
    }
    fpInstance = nullptr;
}

//_____________________________________________________________________________________

G4DNAMolecularReactionTable::G4DNAMolecularReactionTable()
    : G4ITReactionTable()
    , fVerbose(false)
    , fGeometry(nullptr)
    , fpMessenger(new G4ReactionTableMessenger(this))
{
}

G4DNAMolecularReactionTable::~G4DNAMolecularReactionTable() = default;

void G4DNAMolecularReactionTable::SetReaction(G4DNAMolecularReactionData* pReactionData)
{
    const auto pReactant1 = pReactionData->GetReactant1();
    const auto pReactant2 = pReactionData->GetReactant2();

    fReactionData[pReactant1][pReactant2] = pReactionData;
    fReactantsMV[pReactant1].push_back(pReactant2);
    fReactionDataMV[pReactant1].push_back(pReactionData);

    if (pReactant1 != pReactant2)
    {
        fReactionData[pReactant2][pReactant1] = pReactionData;
        fReactantsMV[pReactant2].push_back(pReactant1);
        fReactionDataMV[pReactant2].push_back(pReactionData);
    }

    fVectorOfReactionData.emplace_back(pReactionData);
    pReactionData->SetReactionID((G4int)fVectorOfReactionData.size());
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::SetReaction(G4double reactionRate,
                                              Reactant* pReactant1,
                                              Reactant* pReactant2)
{
    auto reactionData = new G4DNAMolecularReactionData(reactionRate, pReactant1, pReactant2);
    SetReaction(reactionData);
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::PrintTable(G4VDNAReactionModel* pReactionModel)
{
    // Print Reactions and Interaction radius for jump step = 3ps

    G4IosFlagsSaver iosfs(G4cout);

   if (pReactionModel && !(pReactionModel->GetReactionTable()))
   {
       pReactionModel->SetReactionTable(this);
   }

    ReactivesMV::iterator itReactives;

    std::map<Reactant*, std::map<Reactant*, G4bool>> alreadyPrint;

    G4cout << "Number of chemical species involved in reactions = "
           << fReactantsMV.size() << G4endl;

    std::size_t nbPrintable = fReactantsMV.size() * fReactantsMV.size();

    G4String* outputReaction = new G4String[nbPrintable];
    G4String* outputReactionRate = new G4String[nbPrintable];
    G4String* outputRange = new G4String[nbPrintable];
    G4int n = 0;

    for (itReactives = fReactantsMV.begin(); itReactives != fReactantsMV.end();
         ++itReactives)
    {
        Reactant* moleculeA = (Reactant*)itReactives->first;
        const vector<Reactant*>* reactivesVector = CanReactWith(moleculeA);

        if (pReactionModel) pReactionModel->InitialiseToPrint(moleculeA);

        G4int nbReactants = (G4int)fReactantsMV[itReactives->first].size();

        for (G4int iReact = 0; iReact < nbReactants; iReact++)
        {
            auto moleculeB = (Reactant*)(*reactivesVector)[iReact];

            Data* reactionData = fReactionData[moleculeA][moleculeB];

            //-----------------------------------------------------------
            // Name of the reaction
            if (!alreadyPrint[moleculeA][moleculeB])
            {
                outputReaction[n] = moleculeA->GetName() + " + " + moleculeB->GetName();

                G4int nbProducts = reactionData->GetNbProducts();

                if (nbProducts)
                {
                    outputReaction[n] += " -> " + reactionData->GetProduct(0)->GetName();

                    for (G4int j = 1; j < nbProducts; j++)
                    {
                        outputReaction[n] += " + " + reactionData->GetProduct(j)->GetName();
                    }
                }
                else
                {
                    outputReaction[n] += " -> No product";
                }

                //-----------------------------------------------------------
                // Interaction Rate
                outputReactionRate[n] = G4UIcommand::ConvertToString(
                    reactionData->GetObservedReactionRateConstant() / (1e-3 * m3 / (mole * s)));

                //-----------------------------------------------------------
                // Calculation of the Interaction Range
                G4double interactionRange = -1;
                if (pReactionModel) interactionRange =
                    pReactionModel->GetReactionRadius(iReact);

                if (interactionRange != -1)
                {
                    outputRange[n] = G4UIcommand::ConvertToString(
                        interactionRange / nanometer);
                }
                else
                {
                    outputRange[n] = "";
                }

                alreadyPrint[moleculeB][moleculeA] = TRUE;
                n++;
            }
        }
    }
    // G4cout<<"Number of possible reactions: "<< n << G4endl;

    ////////////////////////////////////////////////////////////////////
    // Tableau dynamique en fonction du nombre de caractere maximal dans
    // chaque colonne
    ////////////////////////////////////////////////////////////////////

    G4int maxlengthOutputReaction = -1;
    G4int maxlengthOutputReactionRate = -1;

    for (G4int i = 0; i < n; ++i)
    {
        if (maxlengthOutputReaction < (G4int)outputReaction[i].length())
        {
            maxlengthOutputReaction = (G4int)outputReaction[i].length();
        }
        if (maxlengthOutputReactionRate < (G4int)outputReactionRate[i].length())
        {
            maxlengthOutputReactionRate = (G4int)outputReactionRate[i].length();
        }
    }

    maxlengthOutputReaction += 2;
    maxlengthOutputReactionRate += 2;

    if (maxlengthOutputReaction < 10) maxlengthOutputReaction = 10;
    if (maxlengthOutputReactionRate < 30) maxlengthOutputReactionRate = 30;

    G4String* title;

    if (pReactionModel) title = new G4String[3];
    else title = new G4String[2];

    title[0] = "Reaction";
    title[1] = "Reaction Rate [dm3/(mol*s)]";

    if (pReactionModel) title[2] =
        "Interaction Range for chosen reaction model [nm]";

    G4cout << setfill(' ') << setw(maxlengthOutputReaction) << left << title[0]
        << setw(maxlengthOutputReactionRate) << left << title[1];

    if (pReactionModel) G4cout << setw(2) << left << title[2];

    G4cout << G4endl;

    G4cout.fill('-');
    if (pReactionModel) G4cout.width(
        maxlengthOutputReaction + 2 + maxlengthOutputReactionRate + 2
        + (G4int)title[2].length());
    else G4cout.width(maxlengthOutputReaction + 2 + maxlengthOutputReactionRate);
    G4cout << "-" << G4endl;
    G4cout.fill(' ');

    for (G4int i = 0; i < n; i++)
    {
        G4cout << setw(maxlengthOutputReaction) << left << outputReaction[i]
               << setw(maxlengthOutputReactionRate) << left
               << outputReactionRate[i];

        if (pReactionModel) G4cout << setw(2) << left << outputRange[i];

        G4cout << G4endl;

        G4cout.fill('-');
        if (pReactionModel) G4cout.width(
            maxlengthOutputReaction + 2 + maxlengthOutputReactionRate + 2
            + (G4int)title[2].length());
        else G4cout.width(
            maxlengthOutputReaction + 2 + maxlengthOutputReactionRate);
        G4cout << "-" << G4endl;
        G4cout.fill(' ');
    }

    delete[] title;
    delete[] outputReaction;
    delete[] outputReactionRate;
    delete[] outputRange;
}

//______________________________________________________________________________
// Get/Set methods

G4VDNAMolecularGeometry* G4DNAMolecularReactionTable::GetGeometry() const
{
  return fGeometry;
}

G4DNAMolecularReactionTable::Data*
G4DNAMolecularReactionTable::GetReactionData(Reactant* pReactant1,
                                             Reactant* pReactant2) const
{
    if (fReactionData.empty())
    {
        G4String errMsg = "No reaction table was implemented";
        G4Exception("G4MolecularInteractionTable::GetReactionData", "",
                    FatalErrorInArgument, errMsg);
    }

    auto it1 = fReactionData.find(pReactant1);

    if (it1 == fReactionData.end())
    {
        G4String errMsg =
            "No reaction table was implemented for this molecule Definition : " + pReactant1
            ->GetName();
        G4Exception("G4MolecularInteractionTable::GetReactionData", "",
                    FatalErrorInArgument, errMsg);
    }

    ReactionDataMap::mapped_type::const_iterator it2 = it1->second.find(pReactant2);

    if (it2 == it1->second.end())
    {
        G4cout << "Name : " << pReactant2->GetName() << G4endl;
        G4String errMsg = "No reaction table was implemented for this molecule : "
            + pReactant2->GetName();
        G4Exception("G4MolecularInteractionTable::GetReactionData", "", FatalErrorInArgument, errMsg);
    }

    return (it2->second);
}

const G4DNAMolecularReactionTable::ReactionDataMap& G4DNAMolecularReactionTable::GetAllReactionData()
{
    return fReactionData;
}

G4DNAMolecularReactionTable::DataList G4DNAMolecularReactionTable::GetVectorOfReactionData()
{
    DataList dataList;

    for (const auto& pData : fVectorOfReactionData)
    {
        dataList.emplace_back(pData.get());
    }

    return dataList;
}

//______________________________________________________________________________

const G4DNAMolecularReactionTable::ReactantList*
G4DNAMolecularReactionTable::CanReactWith(Reactant* pMolecule) const
{
    if (fReactantsMV.empty())
    {
        G4String errMsg = "No reaction table was implemented";
        G4Exception("G4MolecularInteractionTable::CanReactWith", "",
                    FatalErrorInArgument, errMsg);
        return 0;
    }

    auto itReactivesMap = fReactantsMV.find(pMolecule);

    if (itReactivesMap == fReactantsMV.end())
    {
#ifdef G4VERBOSE
        if (fVerbose)
        {
            G4String errMsg = "No reaction table was implemented for this molecule : "
                + pMolecule->GetName();
            //        	G4Exception("G4MolecularInteractionTable::CanReactWith","",FatalErrorInArgument, errMsg);
            G4cout << "--- G4MolecularInteractionTable::GetReactionData ---" << G4endl;
            G4cout << errMsg << G4endl;
        }
#endif
        return nullptr;
    }

    if (fVerbose)
    {
        G4cout << " G4MolecularInteractionTable::CanReactWith :" << G4endl;
        G4cout << "You are checking reactants for : " << pMolecule->GetName() << G4endl;
        G4cout << " the number of reactants is : " << itReactivesMap->second.size() << G4endl;

        auto itProductsVector = itReactivesMap->second.cbegin();

        for (; itProductsVector != itReactivesMap->second.end(); itProductsVector++)
        {
            G4cout << (*itProductsVector)->GetName() << G4endl;
        }
    }
    return &(itReactivesMap->second);
}

//______________________________________________________________________________

const G4DNAMolecularReactionTable::SpecificDataList*
G4DNAMolecularReactionTable::GetReativesNData(const G4MolecularConfiguration* molecule) const
{
    if (fReactionData.empty())
    {
        G4String errMsg = "No reaction table was implemented";
        G4Exception("G4MolecularInteractionTable::CanInteractWith", "",
                    FatalErrorInArgument, errMsg);
    }

    ReactionDataMap::const_iterator itReactivesMap = fReactionData.find(molecule);

    if (itReactivesMap == fReactionData.end())
    {
        return nullptr;
    }

    if (fVerbose)
    {
        G4cout << " G4MolecularInteractionTable::CanReactWith :" << G4endl;
        G4cout << "You are checking reactants for : " << molecule->GetName() << G4endl;
        G4cout << " the number of reactants is : " << itReactivesMap->second.size() << G4endl;

        SpecificDataList::const_iterator itProductsVector = itReactivesMap->second.begin();

        for (; itProductsVector != itReactivesMap->second.end(); itProductsVector++)
        {
            G4cout << itProductsVector->first->GetName() << G4endl;
        }
    }
    return &(itReactivesMap->second);
}

//______________________________________________________________________________

const G4DNAMolecularReactionTable::DataList*
G4DNAMolecularReactionTable::GetReactionData(const G4MolecularConfiguration* molecule) const
{
    if (fReactionDataMV.empty())
    {
        G4String errMsg = "No reaction table was implemented";
        G4Exception("G4MolecularInteractionTable::CanInteractWith", "",
                    FatalErrorInArgument, errMsg);
    }
    auto it = fReactionDataMV.find(molecule);

    if (it == fReactionDataMV.end())
    {
        G4String errMsg = "No reaction table was implemented for this molecule Definition : "
            + molecule->GetName();
        G4Exception("G4MolecularInteractionTable::GetReactionData", "", FatalErrorInArgument, errMsg);
    }

    return &(it->second);
}

//______________________________________________________________________________

G4DNAMolecularReactionTable::Data* G4DNAMolecularReactionTable::GetReactionData(const G4String& mol1,
                                                                                const G4String& mol2) const
{
    const auto pConf1 = G4MoleculeTable::GetMoleculeTable()->GetConfiguration(mol1);
    const auto pConf2 = G4MoleculeTable::GetMoleculeTable()->GetConfiguration(mol2);
    return GetReactionData(pConf1, pConf2);
}

//______________________________________________________________________________

void
G4DNAMolecularReactionData::SetPolynomialParameterization(const std::vector<double>& P)
{
    fRateParam = std::bind(PolynomialParam, std::placeholders::_1, P);
}

//______________________________________________________________________________

void G4DNAMolecularReactionData::SetArrehniusParameterization(double A0,
                                                              double E_R)
{
    std::vector<double> P = { A0, E_R };
    fRateParam = std::bind(ArrehniusParam, std::placeholders::_1, P);
}

//______________________________________________________________________________

void G4DNAMolecularReactionData::SetScaledParameterization(double temperature_K,
                                                           double rateCste)
{
    fRateParam = std::bind(ScaledParameterization,
                           std::placeholders::_1,
                           temperature_K,
                           rateCste);
}

//______________________________________________________________________________

void G4DNAMolecularReactionTable::ScaleReactionRateForNewTemperature(double temp_K)
{
    for (const auto& pData : fVectorOfReactionData)
    {
        const_cast<G4DNAMolecularReactionData*>(pData.get())->ScaleForNewTemperature(temp_K);
    }
}

//______________________________________________________________________________

void G4DNAMolecularReactionData::ScaleForNewTemperature(double temp_K)
{
    if (fRateParam)
    {
        SetObservedReactionRateConstant(fRateParam(temp_K));
    }
}

//______________________________________________________________________________

G4DNAMolecularReactionTable::Data*
G4DNAMolecularReactionTable::GetReaction(int reactionID) const
{
    for (auto& pData : fVectorOfReactionData)
    {
        if (pData->GetReactionID() == reactionID)
        {
            return pData.get();
        }
    }
    return nullptr;
}

size_t G4DNAMolecularReactionTable::GetNReactions() const
{
    return fVectorOfReactionData.size();
}

void G4DNAMolecularReactionTable::Reset()
{
    fReactionData.clear();
    fReactantsMV.clear();
    fReactionDataMV.clear();
    fVectorOfReactionData.clear();
}
