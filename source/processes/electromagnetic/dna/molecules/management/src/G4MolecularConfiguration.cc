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
#include "G4MoleculeTable.hh"
#include "G4Serialize.hh"
#include <fstream>

using CLHEP::m2;
using CLHEP::s;
using CLHEP::kelvin;

using namespace std;

#if defined ( WIN32 )
#define __func__ __FUNCTION__
#endif

/*G4ThreadLocal*/G4double G4MolecularConfiguration::fgTemperature = 298; // 310*kelvin;
// 25°C, used to shoot an energy

//______________________________________________________________________________
// G4MolecularConfigurationManager
typedef G4MolecularConfiguration::G4MolecularConfigurationManager
    MolecularConfigurationManager;

MolecularConfigurationManager* G4MolecularConfiguration::fgManager = 0;

G4Mutex MolecularConfigurationManager::fManagerCreationMutex;

int G4MolecularConfiguration::GetNumberOfSpecies()
{
  return GetManager()->GetNumberOfCreatedSpecies();
}

double G4MolecularConfiguration::ReturnDefaultDiffCoeff(const G4Material*,
                                     double,
                                     const G4MolecularConfiguration*
                                     molConf)
{
  return molConf->fDynDiffusionCoefficient;
}

G4MolecularConfiguration::G4MolecularConfiguration(const G4MoleculeDefinition* moleculeDef,
                                                   const G4String& label,
                                                   int charge)
{
  fMoleculeDefinition = moleculeDef;

  fLabel = new G4String(label);

  fMoleculeID = GetManager()->Insert(moleculeDef,
                                     label,
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

  fDiffParam = &G4MolecularConfiguration::ReturnDefaultDiffCoeff;
  fIsFinalized = false;
}

void G4MolecularConfiguration::MakeExceptionIfFinalized()
{
  if(fIsFinalized)
  {
    G4ExceptionDescription errMsg;
    errMsg << "This molecular configuration " << GetName()
           << " is already finalized. Therefore its "
           " properties cannot be changed.";
    G4Exception("G4MolecularConfiguration::MakeExceptionIfFinalized",
                "CONF_FINALIZED",FatalException,errMsg);
  }
}

//______________________________________________________________________________

G4MolecularConfiguration::G4MolecularConfigurationManager*
G4MolecularConfiguration::GetManager()
{
  if (!fgManager)
  {
    G4AutoLock lock(&MolecularConfigurationManager::fManagerCreationMutex);
    if (!fgManager) // double check for MT
    {
      fgManager = new G4MolecularConfiguration::
          G4MolecularConfigurationManager();
    }
    lock.unlock();
  }

  return fgManager;
}

//______________________________________________________________________________

G4MolecularConfiguration::
G4MolecularConfigurationManager::~G4MolecularConfigurationManager()
{
//  G4cout << "Does G4AllocatorList exists= ";
//  G4cout << (G4AllocatorList::GetAllocatorListIfExist() ? "true":"false")
//      << G4endl;

  G4MolecularConfigurationManager::MolElectronConfTable::iterator it1;
  G4MolecularConfigurationManager::ElectronOccupancyTable::
    iterator it2;

  for (it1 = fElecOccTable.begin(); it1 != fElecOccTable.end(); it1++)
  {
    for (it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
    {
      if (it2->second)
      {
        delete it2->second;
      }
    }
  }
  fElecOccTable.clear();
  fgManager = 0;
}

//______________________________________________________________________________
// G4MolecularConfigurationManager
G4int G4MolecularConfiguration::
G4MolecularConfigurationManager::
Insert(const G4MoleculeDefinition* molDef,
                             const G4ElectronOccupancy& eOcc,
                             G4MolecularConfiguration* molConf)
{
  //G4AutoLock lock(&fMoleculeCreationMutex);

  ElectronOccupancyTable& table2 = fElecOccTable[molDef];
  ElectronOccupancyTable::iterator it = table2.find(eOcc);

  if(it == table2.end())
  {
    table2[eOcc] = molConf;
  }
  else
  {
    G4ExceptionDescription errMsg;
    errMsg << "The same molecular configuration seemed to be recorded twice";
    G4Exception("G4MolecularConfigurationManager::"
                "SetMolecularConfiguration(const G4MoleculeDefinition* molDef,"
                "const G4ElectronOccupancy& eOcc,"
                "G4MolecularConfiguration* molConf)",
                "",
                FatalException,
                errMsg
                );
  }

  fLastMoleculeID++;

  fMolConfPerID.push_back(molConf);

  //lock.unlock();
  return fLastMoleculeID;
}

//______________________________________________________________________________

const G4ElectronOccupancy*
G4MolecularConfiguration::G4MolecularConfigurationManager::
FindCommonElectronOccupancy(const G4MoleculeDefinition* molDef,
                            const G4ElectronOccupancy& eOcc)
{
  //G4AutoLock lock(&fMoleculeCreationMutex);

  MolElectronConfTable::iterator it1 = fElecOccTable.find(molDef);

  if(it1 == fElecOccTable.end())
  {
    // TODO = handle exception ?
    return 0;
  }

  ElectronOccupancyTable& table2 = it1->second;
  ElectronOccupancyTable::iterator it2 = table2.find(eOcc);

  //lock.unlock();

  if (it2 == table2.end())
  {
    // TODO = handle exception ?
    return 0;
  }

  return &(it2->first);
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          const G4ElectronOccupancy& eOcc)
{
  MolElectronConfTable::iterator it1 = fElecOccTable.find(molDef);

  if(it1 == fElecOccTable.end()) return 0;

  ElectronOccupancyTable& table2 = it1->second;
  ElectronOccupancyTable::iterator it = table2.find(eOcc);

  if(it == table2.end())
  {
    return 0;
  }
  else
  {
    return it->second;
  }

  return 0;
}

//______________________________________________________________________________

G4int G4MolecularConfiguration::G4MolecularConfigurationManager::
Insert(const G4MoleculeDefinition* molDef,
       int charge,
       G4MolecularConfiguration* molConf)
{

  //G4AutoLock lock(&fMoleculeCreationMutex);
  ChargeTable& table2 = fChargeTable[molDef];
  ChargeTable::iterator it = table2.find(charge);

  if(it == table2.end())
  {
    table2[charge] = molConf;
  }
  else
  {
    //lock.unlock();
    G4ExceptionDescription errMsg;
    errMsg << "The same molecular configuration seemed to be recorded twice";
    G4Exception("G4MolecularConfigurationManager::"
                "SetMolecularConfiguration(const G4MoleculeDefinition* molDef,"
                "int charge,"
                "G4MolecularConfiguration* molConf)",
                "", FatalException, errMsg);
  }

  fLastMoleculeID++;
  fMolConfPerID.push_back(molConf);
  //lock.unlock();
  return fLastMoleculeID;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          int charge)
{
  //G4AutoLock lock(&fMoleculeCreationMutex);

  MolChargeConfTable::iterator it1 = fChargeTable.find(molDef);

  if(it1 == fChargeTable.end()) return 0;

  ChargeTable& table2 = it1->second;
  ChargeTable::iterator it = table2.find(charge);

  if(it == table2.end())
  {
    return 0;
  }
  else
  {
    return it->second;
  }

  return 0;

  //lock.unlock();
  return 0;
}

//______________________________________________________________________________
// Static method in G4MolecularConfiguration
G4MolecularConfiguration*
G4MolecularConfiguration::
GetOrCreateMolecularConfiguration(const G4MoleculeDefinition* molDef)
{
  if (molDef->GetGroundStateElectronOccupancy())
  {
    const G4ElectronOccupancy& elecOcc =
        *molDef->GetGroundStateElectronOccupancy();
    G4MolecularConfiguration* molConf =
        GetManager()->GetMolecularConfiguration(molDef, elecOcc);

    if (molConf)
    {
      return molConf;
    }
    else
    {
      G4MolecularConfiguration* newConf =
          new G4MolecularConfiguration(molDef,
                                       elecOcc);
      newConf->SetUserID(molDef->GetName());
      return newConf;
    }
  }
  else
  {
    G4MolecularConfiguration* molConf =
        GetManager()->GetMolecularConfiguration(molDef, molDef->GetCharge());
    if(molConf)
    {
      return molConf;
    }
    else
    {
      G4MolecularConfiguration* newConf =
          new G4MolecularConfiguration(molDef, molDef->GetCharge());
      newConf->SetUserID(molDef->GetName());
      return newConf;
    }
  }
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
GetOrCreateMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                  const G4ElectronOccupancy& elecOcc)
{
  return GetManager()->GetOrCreateMolecularConfiguration(molDef, elecOcc);

//  G4MolecularConfiguration* molConf =
//      GetManager()->GetMolecularConfiguration(molDef, elecOcc);
//
//  if (molConf)
//  {
//    return molConf;
//  }
//  else
//  {
//    G4MolecularConfiguration* newConf =
//        new G4MolecularConfiguration(molDef, elecOcc);
//    return newConf;
//  }
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
GetOrCreateMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                  int charge)
{
  G4MolecularConfiguration* molConf =
      GetManager()->GetMolecularConfiguration(molDef, charge);

  if(molConf)
  {
    return molConf;
  }
  else
  {
    G4MolecularConfiguration* newConf =
        new G4MolecularConfiguration(molDef, charge);
    return newConf;
  }
}

//______________________________________________________________________________

void G4MolecularConfiguration::DeleteManager()
{
  G4AutoLock lock(&MolecularConfigurationManager::fManagerCreationMutex);
  if (fgManager) delete fgManager;
  fgManager = 0;
  lock.unlock();
}

//______________________________________________________________________________
// G4MolecularConfiguration
G4MolecularConfiguration::
G4MolecularConfiguration(const G4MoleculeDefinition* moleculeDef,
                         const G4ElectronOccupancy& elecOcc,
                         const G4String& label)
{
  fMoleculeDefinition = moleculeDef;

  fMoleculeID = GetManager()->Insert(moleculeDef,
                                     elecOcc,
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

  fLabel = 0; // let it here

  if(label != "")
  {
    SetLabel(label);
  }

  fDiffParam = &G4MolecularConfiguration::ReturnDefaultDiffCoeff;

  fIsFinalized = false;
}

//______________________________________________________________________________

G4MolecularConfiguration::
G4MolecularConfiguration(const G4MoleculeDefinition* moleculeDef,
                         int charge)
{
  fMoleculeDefinition = moleculeDef;

  fMoleculeID = GetManager()->Insert(moleculeDef,
                                     charge,
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

  fLabel = 0;

  fDiffParam = &G4MolecularConfiguration::ReturnDefaultDiffCoeff;

  fIsFinalized = false;
}

//______________________________________________________________________________

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

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
ChangeConfiguration(const G4ElectronOccupancy& newElectronOccupancy) const
{
  G4MolecularConfiguration* output =
      GetManager()->GetMolecularConfiguration(fMoleculeDefinition,
                                              newElectronOccupancy);

  if (!output)
  {
    output = new G4MolecularConfiguration(fMoleculeDefinition,
                                          newElectronOccupancy);
  }
  return output;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::ChangeConfiguration(int charge) const
{
  G4MolecularConfiguration* output =
      GetManager()->GetMolecularConfiguration(fMoleculeDefinition, charge);

  if (!output)
  {
    output = new G4MolecularConfiguration(fMoleculeDefinition, charge);
  }
  return output;
}

//______________________________________________________________________________

G4MolecularConfiguration&
G4MolecularConfiguration::operator=(G4MolecularConfiguration& /*right*/)
{
//  if (&right == this) return *this;
  return *this;
}

//______________________________________________________________________________

/** Method used in Geant4-DNA to excite water molecules
 */
G4MolecularConfiguration*
G4MolecularConfiguration::ExciteMolecule(G4int ExcitedLevel) const
{
//  MakeExceptionIfFinalized();
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);

  newElectronOccupancy.RemoveElectron(ExcitedLevel, 1);
  newElectronOccupancy.AddElectron(5, 1);

  return ChangeConfiguration(newElectronOccupancy);
}

//______________________________________________________________________________

/** Method used in Geant4-DNA to ionize water molecules
 */
G4MolecularConfiguration*
G4MolecularConfiguration::IonizeMolecule(G4int IonizedLevel) const
{
//  MakeExceptionIfFinalized();
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
    G4Exception("G4MolecularConfiguration::IonizeMolecule",
                "",
                FatalErrorInArgument,
                errMsg);
    PrintState();
  }

  // DEBUG
  // PrintState();

  return ChangeConfiguration(newElectronOccupancy);
}

//______________________________________________________________________________

G4MolecularConfiguration* G4MolecularConfiguration::AddElectron(G4int orbit,
                                                                G4int number) const
{
//  MakeExceptionIfFinalized();
  CheckElectronOccupancy(__func__);
  G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);
  newElectronOccupancy.AddElectron(orbit, number);
  return ChangeConfiguration(newElectronOccupancy);
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::RemoveElectron(G4int orbit,
                                         G4int number) const
{
//  MakeExceptionIfFinalized();
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
    G4Exception("G4MolecularConfiguration::RemoveElectron",
                "",
                JustWarning,
                errMsg);
    PrintState();
  }

  return ChangeConfiguration(newElectronOccupancy);
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::MoveOneElectron(G4int orbitToFree,
                                          G4int orbitToFill) const
{
//  MakeExceptionIfFinalized();
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
    G4Exception("G4MolecularConfiguration::MoveOneElectron",
                "",
                FatalErrorInArgument,
                errMsg);
    PrintState();
  }

  return ChangeConfiguration(newElectronOccupancy);
}

//______________________________________________________________________________

const G4String& G4MolecularConfiguration::GetName() const
{
  return fName;
}

//______________________________________________________________________________

const G4String& G4MolecularConfiguration::GetFormatedName() const
{
  return fFormatedName;
}

//______________________________________________________________________________

G4int G4MolecularConfiguration::GetAtomsNumber() const
{
  return fMoleculeDefinition->GetAtomsNumber();
}

//______________________________________________________________________________

G4double G4MolecularConfiguration::GetNbElectrons() const
{
  CheckElectronOccupancy(__func__);
  return fElectronOccupancy->GetTotalOccupancy();
}

//______________________________________________________________________________

void G4MolecularConfiguration::PrintState() const
{
  G4cout << "-------------- Start Printing State " << GetName()
         << " ---------------" << G4endl;

  if (fElectronOccupancy)
  {
    G4cout << "--------------Print electronic state of " << GetName()
           << "---------------" << G4endl;
    fElectronOccupancy->DumpInfo();
    if(fElectronOccupancy==fMoleculeDefinition->GetGroundStateElectronOccupancy())
    {
      G4cout<<"At ground state"<<G4endl;
    }
  }
  else
  {
    G4cout << "--- No electron occupancy set up ---" << G4endl;
  }

  G4cout << "Charge :"
         << fDynCharge
         << G4endl;

  if(fLabel)
  {
    G4cout << "Label :"
           << GetLabel()
           << G4endl;
  }
  G4cout  << "-------------- End Of State " << GetName()
          << " -----------------------" << G4endl;
}

//______________________________________________________________________________

// added - to be transformed in a "Decay method"
const vector<const G4MolecularDissociationChannel*>*
  G4MolecularConfiguration::GetDissociationChannels() const
{
  // if (fElectronOccupancy == 0) return 0;
  return fMoleculeDefinition->GetDecayChannels(this);
}

//______________________________________________________________________________

G4int G4MolecularConfiguration::GetFakeParticleID() const
{
  if(fMoleculeDefinition) return fMoleculeDefinition->GetPDGEncoding();
  else G4Exception("G4MolecularConfiguration::GetMoleculeID",
                   "",
                   FatalErrorInArgument,
                   "You should first enter a molecule definition");

  return INT_MAX;
}

//______________________________________________________________________________

const char* removePath(const char* path)
{
  const char* pDelimeter = strrchr(path, '\\');
  if (pDelimeter) path = pDelimeter + 1;

  pDelimeter = strrchr(path, '/');
  if (pDelimeter) path = pDelimeter + 1;

  return path;
}

//______________________________________________________________________________

void G4MolecularConfiguration::CheckElectronOccupancy(const char* function) const
{
  if (fElectronOccupancy == 0)
  {
    G4String functionName(function);
    G4ExceptionDescription description;
    description
        << "No G4ElectronOccupancy was defined for molecule definition : "
        << fMoleculeDefinition->GetName()
        << ". The definition was probably defined using the charge state, "
            "rather than electron state.";

    G4Exception(functionName, "", FatalErrorInArgument, description);
  }
}

//______________________________________________________________________________

void G4MolecularConfiguration::G4MolecularConfigurationManager::
RecordNewlyLabeledConfiguration(G4MolecularConfiguration* molConf)
{
  //G4AutoLock lock(&fMoleculeCreationMutex);

  LabelTable& tmpMap = fLabelTable[molConf->fMoleculeDefinition];

  LabelTable::iterator it = tmpMap.find(*molConf->fLabel);

  if(it == tmpMap.end())
  {
    tmpMap[*(molConf->fLabel)] = molConf;
  }
  else
  {
    G4ExceptionDescription errMsg;
    errMsg << "The same molecular configuration seemed to be recorded twice";
    G4Exception("G4MolecularConfigurationManager::"
                "SetMolecularConfiguration(const G4MoleculeDefinition* molDef,"
                "const G4String& label,"
                "G4MolecularConfiguration* molConf)",
                "", FatalException, errMsg);
  }

  //lock.unlock();
}

void G4MolecularConfiguration::G4MolecularConfigurationManager::AddUserID(const G4String& userID,
                                                                          G4MolecularConfiguration* molecule)
{
  UserIDTable::iterator it = fUserIDTable.find(userID);

  if(it == fUserIDTable.end())
  {
    fUserIDTable[userID] = molecule;
  }
  else if(molecule != it->second)
  {
    // TODO improve exception
    // exception
    G4ExceptionDescription description;
    description << "The user identifier " << userID
                << " was already given in another configuration in the table"
                << G4endl;
  G4Exception("G4MolecularConfiguration::G4MolecularConfigurationManager::AddUserID",
                "CONF_ALREADY_RECORDED",
                FatalException,
                description);
  }
}

//______________________________________________________________________________

void G4MolecularConfiguration::G4MolecularConfigurationManager::
RemoveMolecularConfigurationFromTable(G4MolecularConfiguration* configuration)
{
  MolElectronConfTable::iterator it1 =
      fElecOccTable.find(configuration->GetDefinition());
  MolElectronConfTable::iterator end = fElecOccTable.end();

  if (it1 == end) return;

  std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::
    iterator it2 =
      it1->second.find(*configuration->GetElectronOccupancy());

  if (it2 == it1->second.end()) return;

  it2->second = 0;
//  it1->second.erase(it2);

  configuration->fElectronOccupancy = 0;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                          const G4String& label)
{
  //G4AutoLock lock(&fMoleculeCreationMutex);

  MolLabelConfTable::iterator it1 = fLabelTable.find(molDef);

  if(it1 == fLabelTable.end()) return 0;

  LabelTable& table2 = it1->second;

  LabelTable::iterator it2 = table2.find(label);

  //lock.unlock();

  if(it2 == table2.end()) return 0;
  return it2->second;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetMolecularConfiguration(int moleculeID)
{
  if(moleculeID > (int) fMolConfPerID.size() ||
     moleculeID < 0) return 0;

  return fMolConfPerID[moleculeID];
}

//______________________________________________________________________________

G4int
G4MolecularConfiguration::G4MolecularConfigurationManager::
Insert(const G4MoleculeDefinition* molDef,
                             const G4String& label,
                             G4MolecularConfiguration* molConf)
{
  G4AutoLock lock(&fMoleculeCreationMutex);
  LabelTable& tmpMap = fLabelTable[molDef];
  LabelTable::iterator it = tmpMap.find(label);

  if(it == tmpMap.end())
  {
    fLastMoleculeID++;
    tmpMap[label] = molConf;
    lock.unlock();
  }
  else
  {
    lock.unlock();
    G4ExceptionDescription errMsg;
    errMsg << "The same molecular configuration seemed to be recorded twice";
    G4Exception("G4MolecularConfigurationManager::"
                "SetMolecularConfiguration(const G4MoleculeDefinition* molDef,"
                "const G4String& label,"
                "G4MolecularConfiguration* molConf)",
                "", FatalException, errMsg);
  }

  fMolConfPerID.push_back(molConf);

  return fLastMoleculeID;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                                    const G4String& label)
{
  return GetManager()->GetMolecularConfiguration(molDef, label);
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::GetMolecularConfiguration(int moleculeID)
{
  return GetManager()->GetMolecularConfiguration(moleculeID);
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::CreateMolecularConfiguration(const G4String& userIdentifier,
                                                       const G4MoleculeDefinition* molDef,
                                                       int charge,
                                                       const G4String& label,
                                                       bool& wasAlreadyCreated)
{
  wasAlreadyCreated = false;
  G4MolecularConfiguration* molConf =
      GetManager()->GetMolecularConfiguration(molDef, charge);

  if (molConf)
  {
    if(molConf->fLabel == 0)
    {
      molConf->SetLabel(label);
      G4ExceptionDescription wMsg ;
      wMsg << "The molecular configuration for the definition named "
             << molDef->GetName()
             << " with charge " << charge << " has already been created "
                 "but with NO label";
      G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                  "DOUBLE_CREATION",
                  JustWarning,
                  wMsg);
    }
    else if(*(molConf->fLabel) == "" )
    {
      molConf->SetLabel(label);
    }
    else if(*(molConf->fLabel) != label)
    {
      G4ExceptionDescription errMsg ;
      errMsg << "The molecular configuration for the definition named "
             << molDef->GetName()
             << " with charge " << charge << " has already been created "
                 "but with a different label :"
             << molConf->GetLabel();
      G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                  "DOUBLE_CREATION",
                  FatalErrorInArgument,
                  errMsg);
      // KILL APP
    }

    if(molConf->fUserIdentifier == "")
    {
      molConf->fUserIdentifier = userIdentifier;

      G4ExceptionDescription wMsg ;
      wMsg << "The molecular configuration for the definition named "
             << molDef->GetName()
             << " with label " << label << " has already been created.";
      G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                  "DOUBLE_CREATION",
                  JustWarning,
                  wMsg);
    }
    else if(molConf->fUserIdentifier != userIdentifier)
    {
      G4ExceptionDescription errMsg ;
      errMsg << "The molecular configuration for the definition named "
             << molDef->GetName()
             << " with label " << label << " has already been created "
                 "BUT with a different user ID :"
             << molConf->fUserIdentifier;
      G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                  "DOUBLE_CREATION",
                  FatalErrorInArgument,
                  errMsg);
      // KILL APP
    }

    wasAlreadyCreated = true;
    return molConf;
  }
  else
  {
    G4MolecularConfiguration* newConf =
        new G4MolecularConfiguration(molDef, label, charge);
    newConf->fUserIdentifier = userIdentifier;

    GetManager()->AddUserID(userIdentifier, newConf);

//    G4MoleculeTable::Instance()->RecordMolecularConfiguration(userIdentifier,
//                                                              newConf);
    return newConf;
  }
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
CreateMolecularConfiguration(const G4String& userIdentifier,
                             const G4MoleculeDefinition* molDef,
                             bool& wasAlreadyCreated)
{
  wasAlreadyCreated = false;
  G4MolecularConfiguration* preRegisteredMolConf =
      GetManager()->GetMolecularConfiguration(userIdentifier);

  if(preRegisteredMolConf)
  {
    if(preRegisteredMolConf->GetDefinition() == molDef)
    {
      wasAlreadyCreated = true;
      return preRegisteredMolConf;
    }
  }

  if(molDef->GetGroundStateElectronOccupancy())
  {
    const G4ElectronOccupancy& elecOcc = *molDef
        ->GetGroundStateElectronOccupancy();
    G4MolecularConfiguration* molConf =
        GetManager()->GetMolecularConfiguration(molDef, elecOcc);

    if(molConf)
    {
      if(molConf->fUserIdentifier == "")
      {
        molConf->fUserIdentifier = userIdentifier;
      }
      else if(molConf->fUserIdentifier != userIdentifier)
      {
        G4ExceptionDescription errMsg;
        errMsg << "A molecular configuration for the definition named "
               << molDef->GetName() << " has already been created "
               "and recorded with a different user ID "
               << molConf->fUserIdentifier;
        G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                    "DOUBLE_CREATION",
                    FatalErrorInArgument,
                    errMsg);
      }
// TODO exception
      G4ExceptionDescription errMsg;
      errMsg << "A molecular configuration for the definition named "
             << molDef->GetName() << " has already been created.";
      G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                  "DOUBLE_CREATION",
                  JustWarning,
                  errMsg);
      wasAlreadyCreated = true;
      return molConf;
    }
    else
    {
      // G4cout << "Create molConf for " << molDef->GetName() << G4endl;
      G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef,
                                                                       elecOcc);
      newConf->fUserIdentifier = userIdentifier;

      GetManager()->AddUserID(userIdentifier, newConf);

//      G4MoleculeTable::Instance()->RecordMolecularConfiguration(userIdentifier,
//                                                                newConf);
      return newConf;
    }
  }
  else
  {
    return CreateMolecularConfiguration(userIdentifier,
                                        molDef,
                                        molDef->GetName(),
                                        molDef->GetCharge(),
                                        wasAlreadyCreated);
  }
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
CreateMolecularConfiguration(const G4String& userIdentifier,
                             const G4MoleculeDefinition* molDef,
                             const G4String& label,
                             bool& wasAlreadyCreated)
{
  assert(label != "");
  wasAlreadyCreated = false;

  G4MolecularConfiguration* molConf =
      GetManager()->GetMolecularConfiguration(molDef, label);
  if(molConf)
  {
    if(molConf->fLabel
       && *molConf->fLabel == label)
    {
      wasAlreadyCreated = true;
      return molConf;
    }
    else if(molConf->fLabel == 0)
    {
      wasAlreadyCreated = true;
      molConf->SetLabel(label);
      return molConf;
    }
    else if(*molConf->fLabel == "")
    {
      wasAlreadyCreated = true;
      molConf->SetLabel(label);
      return molConf;
    }

    molConf->PrintState();
    G4ExceptionDescription errMsg ;
    errMsg << "A molecular configuration for the definition named "
           << molDef->GetName()
           << " has already been created "
              "with user ID "
           << molConf->fUserIdentifier << " and label "
           << molConf->GetLabel();
    G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                "DOUBLE_CREATION",
                FatalErrorInArgument,
                errMsg);
    // KILL APP
  }
  else
  {
    G4MolecularConfiguration* newConf =
      new G4MolecularConfiguration(molDef,
                                   label,
                                   molDef->GetCharge());
    newConf->fUserIdentifier = userIdentifier;

    GetManager()->AddUserID(userIdentifier, newConf);

//    G4MoleculeTable::Instance()->
//        RecordMolecularConfiguration(userIdentifier, newConf);
    return newConf;
  }
  return molConf;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
CreateMolecularConfiguration(const G4String& userIdentifier,
                             const G4MoleculeDefinition* molDef,
                             const G4String& label,
                             const G4ElectronOccupancy& eOcc,
                             bool& wasAlreadyCreated)
{
  assert(label != "");
  wasAlreadyCreated = false;

  G4MolecularConfiguration* molConf =
      GetManager()->GetMolecularConfiguration(molDef, eOcc);

  if(molConf)
  {
    if(molConf->GetElectronOccupancy())
    {
      if(*molConf->GetElectronOccupancy() == eOcc)
      {
        if(molConf->fLabel && *molConf->fLabel == label)
        {
          wasAlreadyCreated = true;
          return molConf;
        }
        else if(molConf->fLabel == 0)
        {
          wasAlreadyCreated = true;
          molConf->SetLabel(label);
          return molConf;
        }
        else if(*molConf->fLabel == "")
        {
          wasAlreadyCreated = true;
          molConf->SetLabel(label);
          return molConf;
        }
      }
    }


    molConf->PrintState();
    G4ExceptionDescription errMsg ;
    errMsg << "A molecular configuration for the definition named "
           << molDef->GetName()
           << " has already been created "
              "with user ID "
           << molConf->fUserIdentifier
           << " and possible different electronic state";
    G4Exception("G4MolecularConfiguration::CreateMolecularConfiguration",
                "DOUBLE_CREATION",
                FatalErrorInArgument,
                errMsg);
  }
  else
  {
    G4MolecularConfiguration* newConf =
      new G4MolecularConfiguration(molDef,
                                   eOcc,
                                   label);
    newConf->fUserIdentifier = userIdentifier;

    GetManager()->AddUserID(userIdentifier, newConf);

//    G4MoleculeTable::Instance()->
//        RecordMolecularConfiguration(userIdentifier, newConf);
    return newConf;
  }
  return molConf;
}


//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetOrCreateMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                  const G4ElectronOccupancy& eOcc)
{
  MolElectronConfTable::iterator it1 = fElecOccTable.find(molDef);

  if(it1 == fElecOccTable.end())
  {
    return new G4MolecularConfiguration(molDef, eOcc);
  }

  ElectronOccupancyTable& table2 = it1->second;
  ElectronOccupancyTable::iterator it = table2.find(eOcc);

  if(it == table2.end())
  {
    G4MolecularConfiguration* molConf =
        new G4MolecularConfiguration(molDef, eOcc);
//    molConf->Finalize();
    return molConf;
  }
  else
  {
    return it->second;
  }

  return 0;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::G4MolecularConfigurationManager::
GetOrCreateMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                  int charge)
{
  MolChargeConfTable::iterator it1 = fChargeTable.find(molDef);

  if(it1 == fChargeTable.end())
  {
    G4AutoLock lock(&fMoleculeCreationMutex);

    G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef, charge);
    return newConf ;
  }

  ChargeTable& table2 = it1->second;
  ChargeTable::iterator it = table2.find(charge);

  if(it == table2.end())
  {
    G4AutoLock lock(&fMoleculeCreationMutex);

    G4MolecularConfiguration* newConf =
        new G4MolecularConfiguration(molDef, charge);
//    newConf->Finalize();
    return newConf ;
  }
  else
  {
    return it->second;
  }

  return 0;
}

//______________________________________________________________________________

void G4MolecularConfiguration::Serialize(std::ostream& out)
{
  G4String moleculeName = fMoleculeDefinition->GetName();
  WRITE(out, moleculeName);

//  if(fLabel)
//   out << fLabel;
//  else
//    out << "";
  WRITE(out,fDynDiffusionCoefficient);
  WRITE(out,fDynVanDerVaalsRadius);
  WRITE(out,fDynDecayTime);
  WRITE(out,fDynMass);
  WRITE(out,fDynCharge);
  WRITE(out,fMoleculeID);
  WRITE(out,fFormatedName);
  WRITE(out,fName);
  WRITE(out,fIsFinalized);
}

//______________________________________________________________________________

void G4MolecularConfiguration::Unserialize(std::istream& in)
{
  G4String moleculeName;
  READ(in, moleculeName);
  fMoleculeDefinition =
      G4MoleculeTable::Instance()->GetMoleculeDefinition(moleculeName);

//  G4String label;
//
//  in.read((char*)(&label), sizeof(label));
//
//  if(label)
//   fLabel = new G4String(label);
//  else
//    fLabel = 0;
  READ(in,fDynDiffusionCoefficient);
  READ(in,fDynVanDerVaalsRadius);
  READ(in,fDynDecayTime);
  READ(in,fDynMass);
  READ(in,fDynCharge);
  READ(in,fMoleculeID);
  READ(in,fFormatedName);
  READ(in,fName);
  READ(in,fIsFinalized);
}

//______________________________________________________________________________

G4MolecularConfiguration* G4MolecularConfiguration::Load(std::istream& in)
{
  return new G4MolecularConfiguration(in);
}

//______________________________________________________________________________

G4MolecularConfiguration::G4MolecularConfiguration(std::istream& in)
{
  fLabel = 0; // TODO: for now not serialized
  Unserialize(in);
  fMoleculeDefinition = 0;
  fElectronOccupancy = 0;
  if(fElectronOccupancy)
  {
    GetManager()->Insert(fMoleculeDefinition, *fElectronOccupancy, this);
    fElectronOccupancy =
        GetManager()->FindCommonElectronOccupancy(fMoleculeDefinition,
                                                  *fElectronOccupancy);

    if(fLabel)
    {
      GetManager()->RecordNewlyLabeledConfiguration(this);
    }
  }
  else if(fLabel)
  {
    fMoleculeID = GetManager()->Insert(fMoleculeDefinition, *fLabel, this);
  }
  else if(fDynCharge)
  {
    fMoleculeID = GetManager()->Insert(fMoleculeDefinition, fDynCharge, this);
  }
}

//______________________________________________________________________________

void G4MolecularConfiguration::SetUserID(const G4String& userID)
{
  fUserIdentifier = userID;
  GetManager()->AddUserID(userID, this);
//  G4MoleculeTable::Instance()->RecordMolecularConfiguration(userID, this);
}

//______________________________________________________________________________

double G4MolecularConfiguration::DiffCoeffWater(double temperature_K)
{
  return pow(10, 4.311
             - 2.722e3/temperature_K
             + 8.565e5/(temperature_K *temperature_K)
             - 1.181e8/(temperature_K*temperature_K*temperature_K ))*1e-9*m2/s;
}

//______________________________________________________________________________

void
G4MolecularConfiguration::
ScaleAllDiffusionCoefficientsOnWater(double temperature_K)
{
  double D_water_0 = DiffCoeffWater(fgTemperature);
  double D_water_f = DiffCoeffWater(temperature_K);

  G4cout << "Scaling factor = " << D_water_f/D_water_0 << G4endl;

  G4ConfigurationIterator it =
      G4MoleculeTable::Instance()->GetConfigurationIterator();

  while(it())
  {
    G4MolecularConfiguration* conf = it.value();
    double D_0 = conf->GetDiffusionCoefficient() ;
    double D_f = D_water_f * D_0 /D_water_0;
    conf->SetDiffusionCoefficient(D_f);
  };
}

//______________________________________________________________________________

void G4MolecularConfiguration::CreateDefaultDiffCoeffParam()
{
  if(bool(fDiffParam) == false)
  {
    fDiffParam = &G4MolecularConfiguration::ReturnDefaultDiffCoeff;
  }
}

//______________________________________________________________________________

void G4MolecularConfiguration::SetGlobalTemperature(G4double temperature)
{
  ScaleAllDiffusionCoefficientsOnWater(temperature);
  fgTemperature = temperature;
}

//______________________________________________________________________________

G4double G4MolecularConfiguration::GetGlobalTemperature()
{
  return fgTemperature;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::
G4MolecularConfigurationManager::GetMolecularConfiguration(const G4String& userID)
{
  for(auto it : fMolConfPerID)
  {
    if(it->GetUserID() == userID) return it;
  }
  return 0;
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MolecularConfiguration::GetMolecularConfiguration(const G4String& userID)
{
  return GetManager()->GetMolecularConfiguration(userID);
}

//______________________________________________________________________________

void G4MolecularConfiguration::FinalizeAll()
{
  const std::vector<G4MolecularConfiguration*>& species =
      GetManager()->GetAllSpecies();

  for(std::size_t i = 0; i < species.size() ; ++i)
  {
    species[i]->Finalize();
  }

}

void G4MolecularConfiguration::PrintAll() //hoang added
{
  const std::vector<G4MolecularConfiguration*>& species =
    GetManager()->GetAllSpecies();
  G4cout<<G4endl;
  G4cout<<"Molecular Config"<<std::setw(25)<<" | Diffusion Coefficient (m2 / s) "<<std::setw(20)<<" | Radius (nm) "<<G4endl;
  G4cout<<"__________________________________________"
            "___________________________________"<<G4endl;
  for(std::size_t i = 0; i < species.size() ; ++i)
  {
    G4cout<<species[i]->GetName()
           <<std::setw(G4int(30 - species[i]->GetName().length()))
           <<right<<species[i]->GetDiffusionCoefficient() * 1.0e3<<std::setw(30)
           <<species[i]->GetVanDerVaalsRadius()/CLHEP::nm<<G4endl;
    G4cout<<"__________________________________________"
              "___________________________________"<<G4endl;
  }

}
