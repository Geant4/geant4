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
// $Id: G4MolecularConfiguration.cc 90769 2015-06-09 10:33:41Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4MolecularConfiguration.hh"
#include "G4MoleculeDefinition.hh"
#include "G4UIcommand.hh"
#include "G4AllocatorList.hh"
#include "G4AutoLock.hh"

using namespace std;

#if defined ( WIN32 )
#define __func__ __FUNCTION__
#endif

//______________________________________________________________
// G4MolecularConfigurationManager
typedef G4MolecularConfiguration::G4MolecularConfigurationManager MolecularConfigurationManager;

MolecularConfigurationManager* G4MolecularConfiguration::fgManager = 0;

G4Mutex MolecularConfigurationManager::fManagerCreationMutex;

G4MolecularConfiguration::G4MolecularConfigurationManager*
G4MolecularConfiguration::GetManager()
{
  if (!fgManager)
  {
    G4AutoLock lock(&MolecularConfigurationManager::fManagerCreationMutex);
    if (!fgManager) // double check for MT
    {
      fgManager = new G4MolecularConfiguration::G4MolecularConfigurationManager;
    }
    lock.unlock();
  }

  return fgManager;
}

G4MolecularConfiguration::
G4MolecularConfigurationManager::~G4MolecularConfigurationManager()
{
//  G4cout << "Does G4AllocatorList exists= ";
//  G4cout << (G4AllocatorList::GetAllocatorListIfExist() ? "true":"false")
//      << G4endl;

  G4MolecularConfigurationManager::MolecularConfigurationTable::iterator it1;
  std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::iterator it2;

  for (it1 = fTable.begin(); it1 != fTable.end(); it1++)
  {
    for (it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
    {
      if (it2->second)
      {
        delete it2->second;
      }
    }
  }
  fTable.clear();
  fgManager = 0;
}

//______________________________________________________________
// G4MolecularConfigurationManager
G4int G4MolecularConfiguration::
G4MolecularConfigurationManager::
SetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          const G4ElectronOccupancy& eOcc,
                          G4MolecularConfiguration* molConf)
{
  G4AutoLock lock(&fMoleculeCreationMutex);
  fTable[molDef][eOcc] = molConf;
  fLastMoleculeID++;
  lock.unlock();
  return fLastMoleculeID;
}

const G4ElectronOccupancy*
G4MolecularConfiguration::G4MolecularConfigurationManager::
FindCommonElectronOccupancy(
    const G4MoleculeDefinition* molDef,
    const G4ElectronOccupancy& eOcc)
{
  std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::iterator it;
  G4AutoLock lock(&fMoleculeCreationMutex);
  it = fTable[molDef].find(eOcc);
  lock.unlock();

  if (it == fTable[molDef].end())
  {
    // TODO = handle exception ?
    return 0;
  }

  return &(it->first);
}

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          const G4ElectronOccupancy& eOcc)
{
  G4AutoLock lock(&fMoleculeCreationMutex);
  G4MolecularConfiguration* output = fTable[molDef][eOcc];
  lock.unlock();
  return output;
}

G4int G4MolecularConfiguration::G4MolecularConfigurationManager::
SetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          int charge,
                          G4MolecularConfiguration* molConf)
{
  G4AutoLock lock(&fMoleculeCreationMutex);
  fChargeTable[molDef][charge] = molConf;
  fLastMoleculeID++;
  lock.unlock();
  return fLastMoleculeID;
}

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          int charge)
{
  G4AutoLock lock(&fMoleculeCreationMutex);
  G4MolecularConfiguration* output = fChargeTable[molDef][charge];
  lock.unlock();
  return output;
}

//______________________________________________________________
// Static method in G4MolecularConfiguration
G4MolecularConfiguration* G4MolecularConfiguration::GetMolecularConfiguration(const G4MoleculeDefinition* molDef)
{
  if (molDef->GetGroundStateElectronOccupancy())
  {
    const G4ElectronOccupancy& elecOcc = *molDef
        ->GetGroundStateElectronOccupancy();
    G4MolecularConfiguration* molConf = GetManager()->GetMolecularConfiguration(
        molDef, elecOcc);

    if (molConf)
    {
      return molConf;
    }
    else
    {
      G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef,
                                                                       elecOcc);
      return newConf;
    }
  }
  else
  {
    return GetMolecularConfiguration(molDef, molDef->GetCharge());
  }
}

G4MolecularConfiguration*
G4MolecularConfiguration::
GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          const G4ElectronOccupancy& elecOcc)
{
  G4MolecularConfiguration* molConf = GetManager()->GetMolecularConfiguration(
      molDef, elecOcc);

  if (molConf)
  {
    return molConf;
  }
  else
  {
    G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef,
                                                                     elecOcc);
    return newConf;
  }
}

G4MolecularConfiguration*
G4MolecularConfiguration::GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                                    int charge)
{
  G4MolecularConfiguration* molConf =
      GetManager()->GetMolecularConfiguration(molDef, charge);

  if (molConf)
  {
    return molConf;
  }
  else
  {
    G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef,
                                                                     charge);
    return newConf;
  }
}

void G4MolecularConfiguration::DeleteManager()
{
  G4AutoLock lock(&MolecularConfigurationManager::fManagerCreationMutex);
  if (fgManager) delete fgManager;
  fgManager = 0;
  lock.unlock();
}

//______________________________________________________________
// G4MolecularConfiguration
G4MolecularConfiguration::G4MolecularConfiguration(const G4MoleculeDefinition* moleculeDef,
                                                   const G4ElectronOccupancy& elecOcc)
{
  fMoleculeDefinition = moleculeDef;

  fMoleculeID = GetManager()->SetMolecularConfiguration(moleculeDef, elecOcc,
                                                        this);
  fElectronOccupancy = GetManager()->FindCommonElectronOccupancy(moleculeDef,
                                                                 elecOcc);

  /*
   fgManager->fTable[fMoleculeDefinition][elecOcc] = this;
   std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::iterator it ;
   it = fgManager->fTable[moleculeDef].find(elecOcc);
   fElectronOccupancy = &(it->first);
   */

  fDynCharge = fMoleculeDefinition->GetNbElectrons()
      - fElectronOccupancy->GetTotalOccupancy()
               + moleculeDef->GetCharge();
  fDynMass = fMoleculeDefinition->GetMass();

  fDynDiffusionCoefficient = fMoleculeDefinition->GetDiffusionCoefficient();
  fDynVanDerVaalsRadius = fMoleculeDefinition->GetVanDerVaalsRadius();
  fDynDecayTime = fMoleculeDefinition->GetDecayTime();

  fName = fMoleculeDefinition->GetName();
  fName += "^";
  fName += G4UIcommand::ConvertToString(fDynCharge);

  fFormatedName = fMoleculeDefinition->GetFormatedName();
  fFormatedName += "^";
  fFormatedName += "{";
  fFormatedName += G4UIcommand::ConvertToString(fDynCharge);
  fFormatedName += "}";
}

G4MolecularConfiguration::G4MolecularConfiguration(const G4MoleculeDefinition* moleculeDef,
                                                   int charge)
{
  fMoleculeDefinition = moleculeDef;

  fMoleculeID = GetManager()->SetMolecularConfiguration(moleculeDef, charge,
                                                        this);
  fElectronOccupancy = 0;

  fDynCharge = charge;
  fDynMass = fMoleculeDefinition->GetMass();

  fDynDiffusionCoefficient = fMoleculeDefinition->GetDiffusionCoefficient();
  fDynVanDerVaalsRadius = fMoleculeDefinition->GetVanDerVaalsRadius();
  fDynDecayTime = fMoleculeDefinition->GetDecayTime();

  fName = fMoleculeDefinition->GetName();
  fName += "^";
  fName += G4UIcommand::ConvertToString(fDynCharge);

  fFormatedName = fMoleculeDefinition->GetFormatedName();
  fFormatedName += "^";
  fFormatedName += "{";
  fFormatedName += G4UIcommand::ConvertToString(fDynCharge);
  fFormatedName += "}";
}

G4MolecularConfiguration::~G4MolecularConfiguration()
{
  if (fgManager) fgManager->RemoveMolecularConfigurationFromTable(this);

//  if (G4AllocatorList::GetAllocatorListIfExist())
//  {
//    if (fElectronOccupancy)
//    {
//      delete fElectronOccupancy;
//      fElectronOccupancy = 0;
//    }
//  }
}

G4MolecularConfiguration* G4MolecularConfiguration::ChangeConfiguration(const G4ElectronOccupancy& newElectronOccupancy)
{
  G4MolecularConfiguration* output = GetManager()->GetMolecularConfiguration(
      fMoleculeDefinition, newElectronOccupancy);

  if (!output)
  {
    output = new G4MolecularConfiguration(fMoleculeDefinition,
                                          newElectronOccupancy);
  }
  return output;
}

G4MolecularConfiguration* G4MolecularConfiguration::ChangeConfiguration(int charge)
{
  G4MolecularConfiguration* output = GetManager()->GetMolecularConfiguration(
      fMoleculeDefinition, charge);

  if (!output)
  {
    output = new G4MolecularConfiguration(fMoleculeDefinition, charge);
  }
  return output;
}

G4MolecularConfiguration& G4MolecularConfiguration::operator=(G4MolecularConfiguration& right)
{
  if (&right == this) return *this;
  return *this;
}

/** Method used in Geant4-DNA to excite water molecules
 */
G4MolecularConfiguration* G4MolecularConfiguration::ExciteMolecule(G4int ExcitedLevel)
{
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);

  newElectronOccupancy.RemoveElectron(ExcitedLevel, 1);
  newElectronOccupancy.AddElectron(5, 1);

  return ChangeConfiguration(newElectronOccupancy);
}

/** Method used in Geant4-DNA to ionize water molecules
 */
G4MolecularConfiguration* G4MolecularConfiguration::IonizeMolecule(G4int IonizedLevel)
{
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);

  if (newElectronOccupancy.GetOccupancy(IonizedLevel) != 0)
  {
    newElectronOccupancy.RemoveElectron(IonizedLevel, 1);
  }
  else
  {
    G4String errMsg = "There is no electron on the orbit "
        + G4UIcommand::ConvertToString(IonizedLevel)
        + " you want to free. The molecule's name you want to ionized is "
        + GetName();
    G4Exception("G4Molecule::IonizeMolecule", "", FatalErrorInArgument, errMsg);
    PrintState();
  }

  // DEBUG
  // PrintState();

  return ChangeConfiguration(newElectronOccupancy);
}

G4MolecularConfiguration* G4MolecularConfiguration::AddElectron(G4int orbit,
                                                                G4int number)
{
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);
  newElectronOccupancy.AddElectron(orbit, number);
  return ChangeConfiguration(newElectronOccupancy);
}

G4MolecularConfiguration* G4MolecularConfiguration::RemoveElectron(G4int orbit,
                                                                   G4int number)
{
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);

  if (newElectronOccupancy.GetOccupancy(orbit) != 0)
  {
    newElectronOccupancy.RemoveElectron(orbit, number);
  }
  else
  {
    G4String errMsg = "There is already no electron into the orbit "
        + G4UIcommand::ConvertToString(orbit)
        + " you want to free. The molecule's name is " + GetName();
    G4Exception("G4Molecule::RemoveElectron", "", JustWarning, errMsg);
    PrintState();
  }

  return ChangeConfiguration(newElectronOccupancy);
}

G4MolecularConfiguration* G4MolecularConfiguration::MoveOneElectron(G4int orbitToFree,
                                                                    G4int orbitToFill)
{
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);

  if (newElectronOccupancy.GetOccupancy(orbitToFree) >= 1)
  {
    newElectronOccupancy.RemoveElectron(orbitToFree, 1);
    newElectronOccupancy.AddElectron(orbitToFill, 1);
  }
  else
  {
    G4String errMsg = "There is no electron on the orbit "
        + G4UIcommand::ConvertToString(orbitToFree)
        + " you want to free. The molecule's name is " + GetName();
    G4Exception("G4Molecule::MoveOneElectron", "", FatalErrorInArgument,
                errMsg);
    PrintState();
  }

  return ChangeConfiguration(newElectronOccupancy);
}

const G4String& G4MolecularConfiguration::GetName() const
{
//  if (fName.isNull())
//  {
//    fName = fMoleculeDefinition->GetName();
//    fName += "^";
//    // fName+= "{";
//    fName += G4UIcommand::ConvertToString(fDynCharge);
//    // fName+= "}";
//  }
  return fName;
}

const G4String& G4MolecularConfiguration::GetFormatedName() const
{
//  if (fFormatedName.isNull())
//  {
//    fFormatedName = fMoleculeDefinition->GetFormatedName();
//    fFormatedName += "^";
//    fFormatedName += "{";
//    fFormatedName += G4UIcommand::ConvertToString(fDynCharge);
//    fFormatedName += "}";
//  }
  return fFormatedName;
}

G4int G4MolecularConfiguration::GetAtomsNumber() const
{
  return fMoleculeDefinition->GetAtomsNumber();
}

G4double G4MolecularConfiguration::GetNbElectrons() const
{
  CheckElectronOccupancy(__func__);
  return fElectronOccupancy->GetTotalOccupancy();
}

void G4MolecularConfiguration::PrintState() const
{
  if (fElectronOccupancy)
  {
    G4cout << "--------------Print electronic state of " << GetName()
           << "---------------" << G4endl;
    fElectronOccupancy->DumpInfo();
    if(fElectronOccupancy==fMoleculeDefinition->GetGroundStateElectronOccupancy())
    {
      G4cout<<"At ground state"<<G4endl;
    }
    else
    {
      if(fMoleculeDefinition->GetDecayTable())
      G4cout<<"Transition :"<<(fMoleculeDefinition->GetDecayTable())->GetExcitedState(fElectronOccupancy)<<G4endl;
    }
  }
  else
  {
    G4cout<<"--- No electron occupancy set up ---"<<G4endl;
  }
}

// added - to be transformed in a "Decay method"
const vector<const G4MolecularDissociationChannel*>* G4MolecularConfiguration::GetDecayChannel() const
{
  if (fElectronOccupancy == 0) return 0;
  return fMoleculeDefinition->GetDecayChannels(fElectronOccupancy);
}

G4int G4MolecularConfiguration::GetFakeParticleID() const
{
  if (fMoleculeDefinition) return fMoleculeDefinition->GetPDGEncoding();
  else G4Exception("G4Molecule::GetMoleculeID", "", FatalErrorInArgument,
                   "You should first enter a molecule defintion");

  return INT_MAX;
}

const char* removePath(const char* path)
{
  const char* pDelimeter = strrchr(path, '\\');
  if (pDelimeter) path = pDelimeter + 1;

  pDelimeter = strrchr(path, '/');
  if (pDelimeter) path = pDelimeter + 1;

  return path;
}

void G4MolecularConfiguration::CheckElectronOccupancy(const char* function) const
{
  if (fElectronOccupancy == 0)
  {
    G4String functionName(function);
    G4ExceptionDescription description;
    description
        << "No G4ElectronOccupancy was defined for molecule definition : "
        << fMoleculeDefinition->GetName()
        << ". The definition was probably defined using the charge state, rather than electron state.";

    G4Exception(functionName, "", FatalErrorInArgument, description);
  }
}

void G4MolecularConfiguration::G4MolecularConfigurationManager::
RemoveMolecularConfigurationFromTable(G4MolecularConfiguration* configuration)
{
  MolecularConfigurationTable::iterator it1 = fTable.find(
      configuration->GetDefinition());
  MolecularConfigurationTable::iterator end = fTable.end();

  if (it1 == end) return;

  std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::iterator it2 =
      it1->second.find(*configuration->GetElectronOccupancy());

  if (it2 == it1->second.end()) return;

  it2->second = 0;
//  it1->second.erase(it2);

  configuration->fElectronOccupancy = 0;
}
