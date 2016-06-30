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
// $Id: G4DNAMolecularReactionTable.cc 95948 2016-03-03 10:40:33Z gcosmo $
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

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::fInstance(0);

G4DNAMolecularReactionData::G4DNAMolecularReactionData() :
    fReactant1(),
    fReactant2(),
    fObservedReactionRate(0.),
    fEffectiveReactionRadius(0.),
    fProducts(0)
{
  fReactionID = 0;
}

//______________________________________________________________________________

G4DNAMolecularReactionData::
G4DNAMolecularReactionData(G4double reactionRate,
                           G4MolecularConfiguration* reactant1,
                           G4MolecularConfiguration* reactant2) :
    fProducts(0)
{
  fObservedReactionRate = reactionRate;
  SetReactant1(reactant1);
  SetReactant2(reactant2);

  G4double sumDiffCoeff(0.);

  if (reactant1 == reactant2)
  {
    sumDiffCoeff = reactant1->GetDiffusionCoefficient();
    fEffectiveReactionRadius = fObservedReactionRate
        / (4 * pi * sumDiffCoeff * Avogadro);
  }
  else
  {
    sumDiffCoeff = reactant1->GetDiffusionCoefficient()
        + reactant2->GetDiffusionCoefficient();
    fEffectiveReactionRadius = fObservedReactionRate /
        (4 * pi * sumDiffCoeff * Avogadro);
  }
  fReactionID = 0;
}

//______________________________________________________________________________

G4DNAMolecularReactionData::
G4DNAMolecularReactionData(G4double reactionRate,
                           const G4String& reactant1,
                           const G4String& reactant2) :
    fProducts(0)
{
  fObservedReactionRate = reactionRate;
  SetReactant1(reactant1);
  SetReactant2(reactant2);

  G4double sumDiffCoeff(0.);

  if (fReactant1 == fReactant2)
  {
    sumDiffCoeff = fReactant1->GetDiffusionCoefficient();
    fEffectiveReactionRadius = fObservedReactionRate
        / (4 * pi * sumDiffCoeff * Avogadro);
  }
  else
  {
    sumDiffCoeff = fReactant1->GetDiffusionCoefficient()
        + fReactant2->GetDiffusionCoefficient();
    fEffectiveReactionRadius = fObservedReactionRate /
        (4 * pi * sumDiffCoeff * Avogadro);
  }
  fReactionID = 0;
}

G4DNAMolecularReactionData::~G4DNAMolecularReactionData()
{
  if (fProducts)
  {
    fProducts->clear();
    delete fProducts;
    fProducts = 0;
  }
}

void G4DNAMolecularReactionData::SetReactant1(G4MolecularConfiguration* reactive)
{
  fReactant1 = reactive;
}


void G4DNAMolecularReactionData::SetReactant2(G4MolecularConfiguration* reactive)
{
  fReactant2 = reactive;
}

void G4DNAMolecularReactionData::SetReactants(G4MolecularConfiguration* reactant1,
                                             G4MolecularConfiguration* reactant2)
{
  fReactant1 = reactant1;
  fReactant2 = reactant2;
}

void G4DNAMolecularReactionData::AddProduct(G4MolecularConfiguration* molecule)
{
  if (!fProducts) fProducts = new std::vector<G4MolecularConfiguration*>();
  fProducts->push_back(molecule);
}

void G4DNAMolecularReactionData::SetReactant1(const G4String& reactive)
{
  fReactant1 = G4MoleculeTable::Instance()->GetConfiguration(reactive);
}
void G4DNAMolecularReactionData::SetReactant2(const G4String& reactive)
{
  fReactant2 = G4MoleculeTable::Instance()->GetConfiguration(reactive);
}
void G4DNAMolecularReactionData::SetReactants(const G4String& reactant1,
                                             const G4String& reactant2)
{
  fReactant1 = G4MoleculeTable::Instance()->GetConfiguration(reactant1);
  fReactant2 = G4MoleculeTable::Instance()->GetConfiguration(reactant2);
}

void G4DNAMolecularReactionData::AddProduct(const G4String& molecule)
{
  if (!fProducts) fProducts = new std::vector<G4MolecularConfiguration*>();
  fProducts->push_back(G4MoleculeTable::Instance()->GetConfiguration(molecule));
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
   return P[0]*G4Exp(P[1]/temp_K)*
          (1e-3 * CLHEP::m3 / (CLHEP::mole * CLHEP::s));
 }

 double G4DNAMolecularReactionData::ScaledParameterization(double temp_K,
     double temp_init,
     double rateCste_init)
 {
   double D0 = G4MolecularConfiguration::DiffCoeffWater(temp_init);
   double Df = G4MolecularConfiguration::DiffCoeffWater(temp_K);
   return Df*rateCste_init/D0;
 }

//==============================================================================
// REACTION TABLE
//==============================================================================

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::GetReactionTable()
{
  if (!fInstance)
  {
    fInstance = new G4DNAMolecularReactionTable();
  }
  return fInstance;
}

//_____________________________________________________________________________________

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::Instance()
{
  if (!fInstance)
  {
    fInstance = new G4DNAMolecularReactionTable();
  }
  return fInstance;
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::DeleteInstance()
{
  // DEBUG
//        G4cout << "G4MolecularReactionTable::DeleteInstance" << G4endl;
  if (fInstance) delete fInstance;
  fInstance = 0;
}

//_____________________________________________________________________________________

G4DNAMolecularReactionTable::G4DNAMolecularReactionTable() :
    G4ITReactionTable()
// ,   fMoleculeHandleManager(G4MoleculeHandleManager::Instance())
{
//    G4cout << "G4DNAMolecularReactionTable::G4DNAMolecularReactionTable()" << G4endl;
  fVerbose = false;
  fpMessenger = new G4ReactionTableMessenger(this);
  return;
}

//_____________________________________________________________________________________

G4DNAMolecularReactionTable::~G4DNAMolecularReactionTable()
{
  // DEBUG
//    G4cout << "G4MolecularReactionTable::~G4MolecularReactionTable" << G4endl;

  if(fpMessenger) delete fpMessenger;

  ReactionDataMap::iterator it1 = fReactionData.begin();
  std::map<G4MolecularConfiguration*,
            const G4DNAMolecularReactionData*>::iterator it2;

  for(; it1 != fReactionData.end(); it1++)
  {
    for(it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
    {
      const G4DNAMolecularReactionData* reactionData = it2->second;
      if(reactionData)
      {
        G4MolecularConfiguration* reactant1 =
            reactionData->GetReactant1();
        G4MolecularConfiguration* reactant2 =
            reactionData->GetReactant2();

        fReactionData[reactant1][reactant2] = 0;
        fReactionData[reactant2][reactant1] = 0;

        delete reactionData;
      }
    }
  }

  fReactionDataMV.clear();
  fReactionData.clear();
  fReactantsMV.clear();
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::SetReaction(G4DNAMolecularReactionData* reactionData)
{
  G4MolecularConfiguration* reactant1 = reactionData->GetReactant1();
  G4MolecularConfiguration* reactant2 = reactionData->GetReactant2();

  fReactionData[reactant1][reactant2] = reactionData;
  fReactantsMV[reactant1].push_back(reactant2);
  fReactionDataMV[reactant1].push_back(reactionData);

  if (reactant1 != reactant2)
  {
    fReactionData[reactant2][reactant1] = reactionData;
    fReactantsMV[reactant2].push_back(reactant1);
    fReactionDataMV[reactant2].push_back(reactionData);
  }

  fVectorOfReactionData.push_back(reactionData);
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::SetReaction(G4double reactionRate,
                                              G4MolecularConfiguration* reactant1,
                                              G4MolecularConfiguration* reactant2)
{
  G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(
      reactionRate, reactant1, reactant2);
  SetReaction(reactionData);
}

//_____________________________________________________________________________________

void G4DNAMolecularReactionTable::PrintTable(G4VDNAReactionModel* pReactionModel)
{
  // Print Reactions and Interaction radius for jump step = 3ps

  G4IosFlagsSaver iosfs(G4cout);

  if (pReactionModel)
  {
    if (!(pReactionModel->GetReactionTable())) pReactionModel->SetReactionTable(
        this);
  }

  ReactivesMV::iterator itReactives;

  map<G4MolecularConfiguration*, map<G4MolecularConfiguration*, G4bool> > alreadyPrint;

  G4cout << "Number of chemical species involved in reactions = "
         << fReactantsMV.size() << G4endl;

  G4int nbPrintable = fReactantsMV.size() * fReactantsMV.size();

  G4String *outputReaction = new G4String[nbPrintable];
  G4String *outputReactionRate = new G4String[nbPrintable];
  G4String *outputRange = new G4String[nbPrintable];
  G4int n = 0;

  for (itReactives = fReactantsMV.begin(); itReactives != fReactantsMV.end();
      itReactives++)
  {
    G4MolecularConfiguration* moleculeA = (G4MolecularConfiguration*) itReactives->first;
    const vector<G4MolecularConfiguration*>* reactivesVector = CanReactWith(moleculeA);

    if (pReactionModel) pReactionModel->InitialiseToPrint(moleculeA);

    G4int nbReactants = fReactantsMV[itReactives->first].size();

    for (G4int iReact = 0; iReact < nbReactants; iReact++)
    {

      G4MolecularConfiguration* moleculeB = (G4MolecularConfiguration*) (*reactivesVector)[iReact];

      const G4DNAMolecularReactionData* reactionData =
          fReactionData[moleculeA][moleculeB];

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

  for (G4int i = 0; i < n; i++)
  {
    if (maxlengthOutputReaction < (G4int) outputReaction[i].length())
    {
      maxlengthOutputReaction = outputReaction[i].length();
    }
    if (maxlengthOutputReactionRate < (G4int) outputReactionRate[i].length())
    {
      maxlengthOutputReactionRate = outputReactionRate[i].length();
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
      + (G4int) title[2].length());
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
        + (G4int) title[2].length());
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

const G4DNAMolecularReactionData*
G4DNAMolecularReactionTable::GetReactionData(G4MolecularConfiguration* reactant1,
                                             G4MolecularConfiguration* reactant2) const
{
  if (fReactionData.empty())
  {
    G4String errMsg = "No reaction table was implemented";
    G4Exception("G4MolecularInteractionTable::GetReactionData", "",
                FatalErrorInArgument, errMsg);
    return 0;
  }

  ReactionDataMap::const_iterator it1 = fReactionData.find(reactant1);

  if (it1 == fReactionData.end())
  {
    G4String errMsg =
        "No reaction table was implemented for this molecule Definition : " + reactant1
            ->GetName();
//  	G4cout << "--- G4MolecularInteractionTable::GetReactionData ---" << G4endl;
//  	G4cout << errMsg << G4endl;
    G4Exception("G4MolecularInteractionTable::GetReactionData", "",
                FatalErrorInArgument, errMsg);
//  	return 0;
  }

  std::map<G4MolecularConfiguration*,
            const G4DNAMolecularReactionData*>::const_iterator it2 =
                 it1->second.find(reactant2);

  if (it2 == it1->second.end())
  {
    G4cout << "Name : " << reactant2->GetName() << G4endl;
    G4String errMsg = "No reaction table was implemented for this molecule : "
    + reactant2 -> GetName();
    G4Exception("G4MolecularInteractionTable::GetReactionData","",FatalErrorInArgument, errMsg);
  }

  return (it2->second);
}

//______________________________________________________________________________

const std::vector<G4MolecularConfiguration*>*
G4DNAMolecularReactionTable::CanReactWith(G4MolecularConfiguration * aMolecule) const
{
  if (fReactantsMV.empty())
  {
    G4String errMsg = "No reaction table was implemented";
    G4Exception("G4MolecularInteractionTable::CanReactWith", "",
                FatalErrorInArgument, errMsg);
    return 0;
  }

  ReactivesMV::const_iterator itReactivesMap = fReactantsMV.find(aMolecule);

  if (itReactivesMap == fReactantsMV.end())
  {
#ifdef G4VERBOSE
    if (fVerbose)
    {
      G4String errMsg = "No reaction table was implemented for this molecule : "
          + aMolecule->GetName();
//        	G4Exception("G4MolecularInteractionTable::CanReactWith","",FatalErrorInArgument, errMsg);
      G4cout << "--- G4MolecularInteractionTable::GetReactionData ---" << G4endl;
      G4cout << errMsg << G4endl;
    }
#endif
    return 0;
  }
  else
  {
    if(fVerbose)
    {
      G4cout<< " G4MolecularInteractionTable::CanReactWith :"<<G4endl;
      G4cout<<"You are checking reactants for : " << aMolecule->GetName()<<G4endl;
      G4cout<<" the number of reactants is : " << itReactivesMap->second.size()<<G4endl;

      std::vector<G4MolecularConfiguration*>::const_iterator itProductsVector =
      itReactivesMap->second.begin();

      for(; itProductsVector != itReactivesMap->second.end(); itProductsVector++)
      {
        G4cout<<(*itProductsVector)->GetName()<<G4endl;
      }
    }
    return &(itReactivesMap->second);
  }
  return 0;
}

//______________________________________________________________________________

const std::map<G4MolecularConfiguration*, const G4DNAMolecularReactionData*>*
G4DNAMolecularReactionTable::GetReativesNData(G4MolecularConfiguration* molecule) const
{
  if (fReactionData.empty())
  {
    G4String errMsg = "No reaction table was implemented";
    G4Exception("G4MolecularInteractionTable::CanInteractWith", "",
                FatalErrorInArgument, errMsg);
    return 0;
  }

  ReactionDataMap::const_iterator itReactivesMap = fReactionData.find(molecule);

  if (itReactivesMap == fReactionData.end())
  {
    return 0;
//    G4cout << "Nom : " << molecule->GetName() << G4endl;
//    G4String errMsg = "No reaction table was implemented for this molecule Definition : "
//    + molecule -> GetName();
//    G4Exception("G4MolecularInteractionTable::CanReactWith","",FatalErrorInArgument, errMsg);
  }
  else
  {
    if(fVerbose)
    {
      G4cout<< " G4MolecularInteractionTable::CanReactWith :"<<G4endl;
      G4cout<<"You are checking reactants for : " << molecule->GetName()<<G4endl;
      G4cout<<" the number of reactants is : " << itReactivesMap->second.size()<<G4endl;

      std::map<G4MolecularConfiguration*,
      const G4DNAMolecularReactionData*>::const_iterator itProductsVector =
      itReactivesMap->second.begin();

      for(; itProductsVector != itReactivesMap->second.end(); itProductsVector++)
      {
        G4cout<<itProductsVector->first->GetName()<<G4endl;
      }
    }
    return &(itReactivesMap->second);
  }

  return 0;
}

//______________________________________________________________________________

const std::vector<const G4DNAMolecularReactionData*>*
G4DNAMolecularReactionTable::GetReactionData(G4MolecularConfiguration* molecule) const
{
  if (fReactionDataMV.empty())
  {
    G4String errMsg = "No reaction table was implemented";
    G4Exception("G4MolecularInteractionTable::CanInteractWith", "",
                FatalErrorInArgument, errMsg);
    return 0;
  }
  ReactionDataMV::const_iterator it = fReactionDataMV.find(molecule);

  if (it == fReactionDataMV.end())
  {
    G4cout << "Nom : " << molecule->GetName() << G4endl;
    G4String errMsg = "No reaction table was implemented for this molecule Definition : "
    + molecule -> GetName();
    G4Exception("G4MolecularInteractionTable::GetReactionData","",FatalErrorInArgument, errMsg);
    return 0; // coverity
  }

  return &(it->second);
}

//______________________________________________________________________________

const G4DNAMolecularReactionData*
G4DNAMolecularReactionTable::GetReactionData(const G4String& mol1,
                                             const G4String& mol2) const
{
  G4MolecularConfiguration* conf1 = G4MoleculeTable::GetMoleculeTable()
      ->GetConfiguration(mol1);
  G4MolecularConfiguration* conf2 = G4MoleculeTable::GetMoleculeTable()
      ->GetConfiguration(mol2);

  return GetReactionData(conf1, conf2);
}

//______________________________________________________________________________

void
G4DNAMolecularReactionData::
SetPolynomialParameterization(const std::vector<double>& P)
{
  fRateParam = std::bind(PolynomialParam, std::placeholders::_1, P);
}

//______________________________________________________________________________

void G4DNAMolecularReactionData::SetArrehniusParameterization(double A0,
                                                              double E_R)
{
  std::vector<double> P = {A0, E_R};

  G4cout << "ici = " << P[0] << G4endl;
  G4cout << "A0 = " << A0 << G4endl;

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
  size_t end = fVectorOfReactionData.size();

  for(size_t i = 0 ; i < end ; ++i)
  {
    ((G4DNAMolecularReactionData*) fVectorOfReactionData[i])->
        ScaleForNewTemperature(temp_K);
  }
}

//______________________________________________________________________________

void G4DNAMolecularReactionData::ScaleForNewTemperature(double temp_K)
{
  if(fRateParam)
  {
    SetObservedReactionRateConstant(fRateParam(temp_K));

//    G4cout <<"PROD RATE = " << fProductionRate << G4endl;
//
//    if(fProductionRate != DBL_MAX && fProductionRate !=0)
//    {
//      SetPartiallyDiffusionControlledReactionByActivation(fObservedReactionRate,
//                                                          fProductionRate);
//    }
  }
}

//______________________________________________________________________________

const G4DNAMolecularReactionData*
G4DNAMolecularReactionTable::GetReaction(int reactionID) const
{
  for(size_t i = 0 ; i < fVectorOfReactionData.size() ; ++i)
  {
    if(fVectorOfReactionData[i]->GetReactionID() == reactionID)
      return fVectorOfReactionData[i];
  }
  return nullptr;
}
