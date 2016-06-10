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
// $Id: G4DNAMolecularReactionTable.cc 90769 2015-06-09 10:33:41Z gcosmo $
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

using namespace std;

class IosFlagSaver
{
public:
  explicit IosFlagSaver(std::ostream& _ios) :
      ios(_ios), f(_ios.flags())
  {
  }
  ~IosFlagSaver()
  {
    ios.flags(f);
  }

//    IosFlagSaver(const IosFlagSaver &rhs) = delete;
//    IosFlagSaver& operator= (const IosFlagSaver& rhs) = delete;

private:
  std::ostream& ios;
  std::ios::fmtflags f;
};

G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::fInstance(0);
//G4ThreadLocal G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::fInstance(0);

G4DNAMolecularReactionData::G4DNAMolecularReactionData() :
    fReactive1(),
    fReactive2(),
    fReactionRate(0.),
    fReducedReactionRadius(0.),
    fProducts(0)
{
  ;
}

G4DNAMolecularReactionData::G4DNAMolecularReactionData(G4double reactionRate,
                                                       const G4Molecule* reactive1,
                                                       const G4Molecule* reactive2) :
    fProducts(0)
{
  fReactionRate = reactionRate;
  SetReactive1(reactive1);
  SetReactive2(reactive2);

  G4double sumDiffCoeff(0.);

  if (*reactive1 == *reactive2)
  {
    sumDiffCoeff = reactive1->GetDiffusionCoefficient();
    fReducedReactionRadius = fReactionRate
        / (4 * pi * sumDiffCoeff * Avogadro);
  }
  else
  {
    sumDiffCoeff = reactive1->GetDiffusionCoefficient()
        + reactive2->GetDiffusionCoefficient();
    fReducedReactionRadius = fReactionRate / (4 * pi * sumDiffCoeff * Avogadro);
  }
}

G4DNAMolecularReactionData::G4DNAMolecularReactionData(G4double reactionRate,
                                                       const G4String& reactive1,
                                                       const G4String& reactive2) :
    fProducts(0)
{
  fReactionRate = reactionRate;
  SetReactive1(reactive1);
  SetReactive2(reactive2);

  G4double sumDiffCoeff(0.);

  if (*fReactive1 == *fReactive2)
  {
    sumDiffCoeff = fReactive1->GetDiffusionCoefficient();
    fReducedReactionRadius = fReactionRate
        / (4 * pi * sumDiffCoeff * Avogadro);
  }
  else
  {
    sumDiffCoeff = fReactive1->GetDiffusionCoefficient()
        + fReactive2->GetDiffusionCoefficient();
    fReducedReactionRadius = fReactionRate / (4 * pi * sumDiffCoeff * Avogadro);
  }
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

void G4DNAMolecularReactionData::SetReactive1(const G4Molecule* reactive)
{
//    fReactive1 = G4MoleculeHandleManager::Instance()->GetMoleculeHandle(reactive);
  fReactive1 = reactive;
}
void G4DNAMolecularReactionData::SetReactive2(const G4Molecule* reactive)
{
//    fReactive2 = G4MoleculeHandleManager::Instance()->GetMoleculeHandle(reactive);
  fReactive2 = reactive;
}
void G4DNAMolecularReactionData::SetReactive(const G4Molecule* reactive1,
                                             const G4Molecule* reactive2)
{
//    fReactive1 = G4MoleculeHandleManager::Instance()->GetMoleculeHandle(reactive1);
//    fReactive2 = G4MoleculeHandleManager::Instance()->GetMoleculeHandle(reactive2);
  fReactive1 = reactive1;
  fReactive2 = reactive2;
}

void G4DNAMolecularReactionData::AddProduct(const G4Molecule* molecule)
{
//    if(!fProducts) fProducts = new std::vector<G4MoleculeHandle>();
  //    fProducts->push_back(G4MoleculeHandleManager::Instance()->GetMoleculeHandle(molecule));
  if (!fProducts) fProducts = new std::vector<const G4Molecule*>();
  fProducts->push_back(molecule);
}

void G4DNAMolecularReactionData::SetReactive1(const G4String& reactive)
{
  fReactive1 = G4MoleculeTable::Instance()->GetMoleculeModel(reactive);
}
void G4DNAMolecularReactionData::SetReactive2(const G4String& reactive)
{
  fReactive2 = G4MoleculeTable::Instance()->GetMoleculeModel(reactive);
}
void G4DNAMolecularReactionData::SetReactive(const G4String& reactive1,
                                             const G4String& reactive2)
{
  fReactive1 = G4MoleculeTable::Instance()->GetMoleculeModel(reactive1);
  fReactive2 = G4MoleculeTable::Instance()->GetMoleculeModel(reactive2);
}

void G4DNAMolecularReactionData::AddProduct(const G4String& molecule)
{
//    if(!fProducts) fProducts = new std::vector<G4MoleculeHandle>();
  if (!fProducts) fProducts = new std::vector<const G4Molecule*>();
  fProducts->push_back(G4MoleculeTable::Instance()->GetMoleculeModel(molecule));
}
//_____________________________________________________________________________________
G4DNAMolecularReactionTable* G4DNAMolecularReactionTable::GetReactionTable()
{
  if (!fInstance)
  {
    fInstance = new G4DNAMolecularReactionTable();
  }
  return fInstance;
}

void G4DNAMolecularReactionTable::DeleteInstance()
{
  // DEBUG
//        G4cout << "G4MolecularReactionTable::DeleteInstance" << G4endl;
  if (fInstance) delete fInstance;
  fInstance = 0;
}
//_____________________________________________________________________________________
G4DNAMolecularReactionTable::G4DNAMolecularReactionTable() :
    G4ITReactionTable(),
    fMoleculeHandleManager(G4MoleculeHandleManager::Instance())
{
//    G4cout << "G4DNAMolecularReactionTable::G4DNAMolecularReactionTable()" << G4endl;
  fVerbose = false;
  return;
}
//_____________________________________________________________________________________
G4DNAMolecularReactionTable::~G4DNAMolecularReactionTable()
{
  // DEBUG
//    G4cout << "G4MolecularReactionTable::~G4MolecularReactionTable" << G4endl;

  /*
   ReactionDataMap::iterator it1 = fReactionData.begin();

   std::map<const G4Molecule*,
   const G4DNAMolecularReactionData*,
   compMoleculeP>::iterator it2;

   for(;it1!=fReactionData.end();it1++)
   {
   for(it2 = it1->second.begin();it2 != it1->second.end();it2++)
   {
   const G4DNAMolecularReactionData* reactionData = it2->second;
   if(reactionData)
   {
   const G4Molecule* reactive1 = reactionData->GetReactive1();
   const G4Molecule* reactive2 = reactionData->GetReactive2();

   fReactionData[reactive1][reactive2] = 0;
   fReactionData[reactive2][reactive1] = 0;

   delete reactionData;
   }
   }
   }
   */
  fReactionDataMV.clear();
  fReactionData.clear();
  fReactivesMV.clear();
}
//_____________________________________________________________________________________
void G4DNAMolecularReactionTable::SetReaction(G4DNAMolecularReactionData* reactionData)
{
  const G4Molecule* reactive1 = reactionData->GetReactive1();
  const G4Molecule* reactive2 = reactionData->GetReactive2();

  fReactionData[reactive1][reactive2] = reactionData;
  fReactivesMV[reactive1].push_back(reactive2);
  fReactionDataMV[reactive1].push_back(reactionData);

  if (reactive1 != reactive2)
  {
    fReactionData[reactive2][reactive1] = reactionData;
    fReactivesMV[reactive2].push_back(reactive1);
    fReactionDataMV[reactive2].push_back(reactionData);
  }
}
//_____________________________________________________________________________________
void G4DNAMolecularReactionTable::SetReaction(G4double reactionRate,
                                              const G4Molecule* reactive1,
                                              const G4Molecule* reactive2)
{
  G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(
      reactionRate, reactive1, reactive2);
  SetReaction(reactionData);
}
//_____________________________________________________________________________________
void G4DNAMolecularReactionTable::PrintTable(G4VDNAReactionModel* pReactionModel)
{
  // Print Reactions and Interaction radius for jump step = 3ps

  IosFlagSaver iosfs(G4cout);

  if (pReactionModel)
  {
    if (!(pReactionModel->GetReactionTable())) pReactionModel->SetReactionTable(
        this);
  }

  ReactivesMV::iterator itReactives;

  map<G4Molecule*, map<G4Molecule*, G4bool> > alreadyPrint;

  G4cout << "Number of chemical species involved in reactions = "
         << fReactivesMV.size() << G4endl;

  G4int nbPrintable = fReactivesMV.size() * fReactivesMV.size();

  G4String *outputReaction = new G4String[nbPrintable];
  G4String *outputReactionRate = new G4String[nbPrintable];
  G4String *outputRange = new G4String[nbPrintable];
  G4int n = 0;

  for (itReactives = fReactivesMV.begin(); itReactives != fReactivesMV.end();
      itReactives++)
  {
    G4Molecule* moleculeA = (G4Molecule*) itReactives->first;
    const vector<const G4Molecule*>* reactivesVector = CanReactWith(moleculeA);

    if (pReactionModel) pReactionModel->InitialiseToPrint(moleculeA);

    G4int nbReactants = fReactivesMV[itReactives->first].size();

    for (G4int iReact = 0; iReact < nbReactants; iReact++)
    {

      G4Molecule* moleculeB = (G4Molecule*) (*reactivesVector)[iReact];

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
            reactionData->GetReactionRate() / (1e-3 * m3 / (mole * s)));

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
//_____________________________________________________________________________________
// Get/Set methods

const G4DNAMolecularReactionData*
G4DNAMolecularReactionTable::GetReactionData(const G4Molecule* reactive1,
                                             const G4Molecule* reactive2) const
{
  if (fReactionData.empty())
  {
    G4String errMsg = "No reaction table was implemented";
    G4Exception("G4MolecularInteractionTable::GetReactionData", "",
                FatalErrorInArgument, errMsg);
    return 0;
  }

  ReactionDataMap::const_iterator it1 = fReactionData.find(reactive1);

  if (it1 == fReactionData.end())
  {
    G4String errMsg =
        "No reaction table was implemented for this molecule Definition : " + reactive1
            ->GetName();
//  	G4cout << "--- G4MolecularInteractionTable::GetReactionData ---" << G4endl;
//  	G4cout << errMsg << G4endl;
    G4Exception("G4MolecularInteractionTable::GetReactionData", "",
                FatalErrorInArgument, errMsg);
//  	return 0;
  }

  std::map<const G4Molecule*, const G4DNAMolecularReactionData*, compMoleculeP>::const_iterator it2 =
      it1->second.find(reactive2);

  if (it2 == it1->second.end())
  {
    G4cout << "Nom : " << reactive2->GetName() << G4endl;
    G4String errMsg = "No reaction table was implemented for this molecule : "
    + reactive2 -> GetName();
    G4Exception("G4MolecularInteractionTable::GetReactionData","",FatalErrorInArgument, errMsg);
  }

  return (it2->second);
}

const std::vector<const G4Molecule*>*
G4DNAMolecularReactionTable::CanReactWith(const G4Molecule * aMolecule) const
{
  if (fReactivesMV.empty())
  {
    G4String errMsg = "No reaction table was implemented";
    G4Exception("G4MolecularInteractionTable::CanReactWith", "",
                FatalErrorInArgument, errMsg);
    return 0;
  }

  ReactivesMV::const_iterator itReactivesMap = fReactivesMV.find(aMolecule);

  if (itReactivesMap == fReactivesMV.end())
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

      std::vector<const G4Molecule*>::const_iterator itProductsVector =
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

//_____________________________________________________________________________________
const std::map<const G4Molecule*, const G4DNAMolecularReactionData*,
    compMoleculeP>*
G4DNAMolecularReactionTable::GetReativesNData(const G4Molecule* molecule) const
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
    G4cout << "Nom : " << molecule->GetName() << G4endl;
    G4String errMsg = "No reaction table was implemented for this molecule Definition : "
    + molecule -> GetName();
    G4Exception("G4MolecularInteractionTable::CanReactWith","",FatalErrorInArgument, errMsg);
  }
  else
  {
    if(fVerbose)
    {
      G4cout<< " G4MolecularInteractionTable::CanReactWith :"<<G4endl;
      G4cout<<"You are checking reactants for : " << molecule->GetName()<<G4endl;
      G4cout<<" the number of reactants is : " << itReactivesMap->second.size()<<G4endl;

      std::map<const G4Molecule*,
      const G4DNAMolecularReactionData*,
      compMoleculeP>::const_iterator itProductsVector =
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

const std::vector<const G4DNAMolecularReactionData*>*
G4DNAMolecularReactionTable::GetReactionData(const G4Molecule* molecule) const
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
