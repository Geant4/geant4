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
#include <memory>
#include "G4DNAScavengerMaterial.hh"
#include "G4StateManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4PhysicalConstants.hh"
#include "G4MolecularConfiguration.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNABoundingBox.hh"
#include "G4VChemistryWorld.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "G4MoleculeTable.hh"

using namespace std;

//------------------------------------------------------------------------------

G4DNAScavengerMaterial::G4DNAScavengerMaterial(
  G4VChemistryWorld* pChemistryInfo)
  : G4VScavengerMaterial()
  , fpChemistryInfo(pChemistryInfo)
  , fIsInitialized(false)
  , fCounterAgainstTime(false)
  , fVerbose(0)
{
  Initialize();
}

//------------------------------------------------------------------------------

void G4DNAScavengerMaterial::Initialize()
{
  if(fIsInitialized)
  {
    return;
  }

  if(fpChemistryInfo->size() == 0)
  {
    G4cout << "G4DNAScavengerMaterial existed but empty" << G4endl;
  }
  Reset();
  fIsInitialized = true;
}

G4double G4DNAScavengerMaterial::GetNumberMoleculePerVolumeUnitForMaterialConf(
  MolType matConf) const
{
  // no change these molecules
  if(G4MoleculeTable::Instance()->GetConfiguration("H2O") == matConf)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "matConf : "<<matConf->GetName();
    G4Exception("G4DNAScavengerMaterial::GetNumberMoleculePerVolumeUnitForMaterialConf", "G4DNAScavengerMaterial001",
                FatalErrorInArgument, exceptionDescription);
  }

  auto iter = fScavengerTable.find(matConf);
  if(iter == fScavengerTable.end())
  {
    return 0;
  }
  else
  {
    if(iter->second >= 1)
    {
      return (floor)(iter->second);
    }
    else
    {
      return 0;
    }
  }
}

void G4DNAScavengerMaterial::ReduceNumberMoleculePerVolumeUnitForMaterialConf(
  MolType matConf, G4double time)
{
  // no change these molecules
  if(G4MoleculeTable::Instance()->GetConfiguration("H2O") == matConf ||
     G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)") ==
       matConf ||  // suppose that pH is not changed during simu
     G4MoleculeTable::Instance()->GetConfiguration("OHm(B)") == matConf)
  {
    // G4cout<<"moletype : "<<matConf->GetName()<<G4endl;
    // kobs is already counted these molecule concentrations
    return;
  }
  if(!find(matConf))  // matConf must greater than 0
  {
    return;
  }
  fScavengerTable[matConf]--;
  if(fScavengerTable[matConf] < 0)  // projection
  {
    assert(false);
  }

  if(fCounterAgainstTime)
  {
    RemoveAMoleculeAtTime(matConf, time);
  }
}

void G4DNAScavengerMaterial::AddNumberMoleculePerVolumeUnitForMaterialConf(
  MolType matConf, G4double time)
{
  // no change these molecules

  if(G4MoleculeTable::Instance()->GetConfiguration("H2O") == matConf ||
     G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)") ==
       matConf ||  // pH has no change
     G4MoleculeTable::Instance()->GetConfiguration("OHm(B)") == matConf)
  {
    // G4cout<<"moletype : "<<matConf->GetName()<<G4endl;
    // kobs is already counted these molecule concentrations
    return;
  }

  auto it = fScavengerTable.find(matConf);
  if(it == fScavengerTable.end())  // matConf must be in fScavengerTable
  {
    return;
  }
  fScavengerTable[matConf]++;

  if(fCounterAgainstTime)
  {
    AddAMoleculeAtTime(matConf, time);
  }
}

void G4DNAScavengerMaterial::PrintInfo()
{
  auto pConfinedBox = fpChemistryInfo->GetChemistryBoundary();
  auto iter         = fpChemistryInfo->begin();
  G4cout << "**************************************************************"
         << G4endl;
  for(; iter != fpChemistryInfo->end(); iter++)
  {
    auto containedConf = iter->first;
    // auto concentration = iter->second;
    auto concentration =
      fScavengerTable[containedConf] / (Avogadro * pConfinedBox->Volume());
    G4cout << "Scavenger:" << containedConf->GetName() << "  : "
           << concentration / 1.0e-6 /*mm3 to L*/ << " (M)  with : "
           << fScavengerTable[containedConf] << " (molecules)"
           << "in: " << pConfinedBox->Volume() / (um * um * um) << " (um3)"
           << G4endl;
    if(fScavengerTable[containedConf] < 1)
    {
      G4cout << "!!!!!!!!!!!!! this molecule has less one molecule for "
                "considered volume"
             << G4endl;
      // assert(false);
    }
    if(fVerbose != 0)
    {
      Dump();
    }
  }
  G4cout << "**************************************************************"
         << G4endl;
}

void G4DNAScavengerMaterial::Reset()
{
  if(fpChemistryInfo == nullptr)
  {
    return;
  }

  if(fpChemistryInfo->size() == 0)
  {
    return;
  }

  fScavengerTable.clear();
  fCounterMap.clear();
  fpLastSearch.reset(nullptr);

  auto pConfinedBox = fpChemistryInfo->GetChemistryBoundary();
  auto iter         = fpChemistryInfo->begin();
  for(; iter != fpChemistryInfo->end(); iter++)
  {
    auto containedConf = iter->first;
    auto concentration = iter->second;
    fScavengerTable[containedConf] =
      floor(Avogadro * concentration * pConfinedBox->Volume());
    fCounterMap[containedConf][1 * picosecond] =
      floor(Avogadro * concentration * pConfinedBox->Volume());
  }
  if(fVerbose != 0){PrintInfo();}
}

//------------------------------------------------------------------------------

void G4DNAScavengerMaterial::AddAMoleculeAtTime(
  MolType molecule, G4double time, const G4ThreeVector* /*position*/,
  int number)
{
  if(fVerbose != 0)
  {
    G4cout << "G4DNAScavengerMaterial::AddAMoleculeAtTime : "
           << molecule->GetName() << " at time : " << G4BestUnit(time, "Time")
           << G4endl;
  }

  auto counterMap_i = fCounterMap.find(molecule);

  if(counterMap_i == fCounterMap.end())
  {
    fCounterMap[molecule][time] = number;
  }
  else if(counterMap_i->second.empty())
  {
    counterMap_i->second[time] = number;
  }
  else
  {
    auto end = counterMap_i->second.rbegin();

    if(end->first <= time || fabs(end->first - time) <=
                               G4::MoleculeCounter::TimePrecision::fPrecision)
    {
      G4double newValue            = end->second + number;
      counterMap_i->second[time] = newValue;
      if(newValue != (floor)(fScavengerTable[molecule]))  // protection
      {
        G4String errMsg = "You are trying to add wrong molecule ";
        G4Exception("AddAMoleculeAtTime", "",
                    FatalErrorInArgument, errMsg);

      }
    }
  }
}

//------------------------------------------------------------------------------

void G4DNAScavengerMaterial::RemoveAMoleculeAtTime(
  MolType pMolecule, G4double time, const G4ThreeVector* /*position*/,
  int number)
{
  NbMoleculeInTime& nbMolPerTime = fCounterMap[pMolecule];

  if(fVerbose != 0)
  {
    auto it_ = nbMolPerTime.rbegin();
    G4cout << "G4DNAScavengerMaterial::RemoveAMoleculeAtTime : "
           << pMolecule->GetName() << " at time : " << G4BestUnit(time, "Time")

           << " form : " << it_->second << G4endl;
  }

  if(nbMolPerTime.empty())
  {
    Dump();
    G4String errMsg = "You are trying to remove molecule " +
                      pMolecule->GetName() +
                      " from the counter while this kind of molecules has not "
                      "been registered yet";
    G4Exception("G4DNAScavengerMaterial::RemoveAMoleculeAtTime", "",
                FatalErrorInArgument, errMsg);

    return;
  }
  else
  {
    auto it = nbMolPerTime.rbegin();

    if(it == nbMolPerTime.rend())
    {
      it--;

      G4String errMsg = "There was no " + pMolecule->GetName() +
                        " recorded at the time or even before the time asked";
      G4Exception("G4DNAScavengerMaterial::RemoveAMoleculeAtTime", "",
                  FatalErrorInArgument, errMsg);
    }

    G4double finalN = it->second - number;
    if(finalN < 0)
    {
      Dump();

      G4cout << "fScavengerTable : " << pMolecule->GetName() << " : "
             << (fScavengerTable[pMolecule]) << G4endl;

      G4ExceptionDescription errMsg;
      errMsg << "After removal of " << number << " species of "
             << " " << it->second << " " << pMolecule->GetName()
             << " the final number at time " << G4BestUnit(time, "Time")
             << " is less than zero and so not valid." << G4endl;
      G4cout << " Global time is "
             << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
             << ". Previous selected time is " << G4BestUnit(it->first, "Time")
             << G4endl;
      G4Exception("G4DNAScavengerMaterial::RemoveAMoleculeAtTime", "N_INF_0",
                  FatalException, errMsg);
    }
    nbMolPerTime[time] = finalN;

    if(finalN != (floor)(fScavengerTable[pMolecule]))  // protection
    {
      assert(false);
    }
  }
}

void G4DNAScavengerMaterial::Dump()
{
  auto pConfinedBox = fpChemistryInfo->GetChemistryBoundary();
  auto V            = pConfinedBox->Volume();
  for(const auto& it : fCounterMap)
  {
    auto pReactant = it.first;

    G4cout << " --- > For " << pReactant->GetName() << G4endl;

    for(const auto& it2 : it.second)
    {
      G4cout << " " << G4BestUnit(it2.first, "Time") << "    "
             << it2.second / (Avogadro * V * 1.0e-6 /*mm3 to L*/) << G4endl;
    }
  }
}

int64_t G4DNAScavengerMaterial::GetNMoleculesAtTime(MolType molecule, G4double time)
{
  if(!fCounterAgainstTime)
  {
    G4cout << "fCounterAgainstTime == false" << G4endl;
    assert(false);
  }

  G4bool sameTypeOfMolecule = SearchTimeMap(molecule);
  auto output = SearchUpperBoundTime(time, sameTypeOfMolecule);
  if(output < 0)
  {
    G4ExceptionDescription errMsg;
    errMsg << "N molecules not valid < 0 : "<<
                      molecule->GetName() <<" N : "<< output << G4endl;
    G4Exception("G4DNAScavengerMaterial::GetNMoleculesAtTime", "",
                FatalErrorInArgument, errMsg);
  }
  return output;
}

G4bool G4DNAScavengerMaterial::SearchTimeMap(MolType molecule)
{
  if(fpLastSearch == nullptr)
  {
    fpLastSearch = std::make_unique<Search>();
  }
  else
  {
    if(fpLastSearch->fLowerBoundSet &&
       fpLastSearch->fLastMoleculeSearched->first == molecule)
    {
      return true;
    }
  }

  auto mol_it                         = fCounterMap.find(molecule);
  fpLastSearch->fLastMoleculeSearched = mol_it;

  if(mol_it != fCounterMap.end())
  {
    fpLastSearch->fLowerBoundTime =
      fpLastSearch->fLastMoleculeSearched->second.end();
    fpLastSearch->fLowerBoundSet = true;
  }
  else
  {
    fpLastSearch->fLowerBoundSet = false;
  }

  return false;
}

//------------------------------------------------------------------------------

int64_t G4DNAScavengerMaterial::SearchUpperBoundTime(G4double time,
                                                 G4bool sameTypeOfMolecule)
{
  auto mol_it = fpLastSearch->fLastMoleculeSearched;
  if(mol_it == fCounterMap.end())
  {
    return 0;
  }

  NbMoleculeInTime& timeMap = mol_it->second;
  if(timeMap.empty())
  {
    return 0;
  }

  if(sameTypeOfMolecule)
  {
    if(fpLastSearch->fLowerBoundSet &&
       fpLastSearch->fLowerBoundTime != timeMap.end())
    {
      if(fpLastSearch->fLowerBoundTime->first < time)
      {
        auto upperToLast = fpLastSearch->fLowerBoundTime;
        upperToLast++;

        if(upperToLast == timeMap.end())
        {
          return fpLastSearch->fLowerBoundTime->second;
        }

        if(upperToLast->first > time)
        {
          return fpLastSearch->fLowerBoundTime->second;
        }
      }
    }
  }

  auto up_time_it = timeMap.upper_bound(time);

  if(up_time_it == timeMap.end())
  {
    auto last_time = timeMap.rbegin();
    return last_time->second;
  }
  if(up_time_it == timeMap.begin())
  {
    return 0;
  }

  up_time_it--;

  fpLastSearch->fLowerBoundTime = up_time_it;
  fpLastSearch->fLowerBoundSet  = true;

  return fpLastSearch->fLowerBoundTime->second;
}