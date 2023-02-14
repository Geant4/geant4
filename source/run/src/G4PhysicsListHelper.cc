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
// G4PhysicsListHelper implementation
//
// Author: H.Kurashige, 29 April 2011
// --------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "G4PhysicsListHelper.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"
#include "G4ios.hh"

#include "G4CoupledTransportation.hh"
#include "G4RunManagerKernel.hh"
#include "G4ScoringManager.hh"
#include "G4Transportation.hh"

#include "G4ProcessManager.hh"

#include "G4ProcessType.hh"
#include "G4TransportationProcessType.hh"
#include "G4DecayProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4OpProcessSubType.hh"
#include "G4HadronicProcessType.hh"


G4ThreadLocal G4PhysicsListHelper* G4PhysicsListHelper::pPLHelper = nullptr;

// --------------------------------------------------------------------
G4PhysicsListHelper::G4PhysicsListHelper()
{
  // pointer to the particle table
  theParticleTable  = G4ParticleTable::GetParticleTable();
  aParticleIterator = theParticleTable->GetIterator();

  ReadOrdingParameterTable();

#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    DumpOrdingParameterTable();
  }
#endif
}

// --------------------------------------------------------------------
G4PhysicsListHelper::~G4PhysicsListHelper()
{
  if(theTable != nullptr)
  {
    theTable->clear();
    delete theTable;
    theTable    = nullptr;
    sizeOfTable = 0;
  }
}

// --------------------------------------------------------------------
G4PhysicsListHelper* G4PhysicsListHelper::GetPhysicsListHelper()
{
  if(pPLHelper == nullptr)
  {
    static G4ThreadLocalSingleton<G4PhysicsListHelper> inst;
    pPLHelper = inst.Instance();
  }
  return pPLHelper;
}

// --------------------------------------------------------------------
void G4PhysicsListHelper::CheckParticleList() const
{
  G4bool isElectron         = false;
  G4bool isPositron         = false;
  G4bool isGamma            = false;
  G4bool isProton           = false;
  G4bool isGenericIon       = false;
  G4bool isAnyIon           = false;
  G4bool isAnyChargedBaryon = false;
  G4bool isEmProc           = false;

  // loop over all particles in G4ParticleTable
  aParticleIterator->reset();
  while((*aParticleIterator)())
  {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String name                  = particle->GetParticleName();
    // check if any EM process exists
    if(!isEmProc)
    {
      G4ProcessVector* list = particle->GetProcessManager()->GetProcessList();
      for(G4int idx = 0; idx < (G4int)list->size(); ++idx)
      {
        isEmProc = ((*list)[idx])->GetProcessType() == fElectromagnetic;
        if(isEmProc)
          break;
      }
    }

    if(name == "e-")
      isElectron = true;
    else if(name == "e+")
      isPositron = true;
    else if(name == "gamma")
      isGamma = true;
    else if(name == "GenericIon")
      isGenericIon = true;
    else if(name == "proton")
      isProton = true;
    else if(particle->GetParticleType() == "nucleus")
      isAnyIon = true;
    else if(particle->GetParticleType() == "baryon")
    {
      if(particle->GetPDGCharge() != 0.0)
        isAnyChargedBaryon = true;
    }
  }

  if(!isEmProc)
    return;

  // RULE 1
  //  e+, e- and gamma should exist
  //   if one of them exist
  G4bool isEmBasic        = isElectron || isPositron || isGamma;
  G4bool isMissingEmBasic = !isElectron || !isPositron || !isGamma;
  if(isEmBasic && isMissingEmBasic)
  {
    G4String missingName = "";
    if(!isElectron)
      missingName += "e- ";
    if(!isPositron)
      missingName += "e+ ";
    if(!isGamma)
      missingName += "gamma ";

#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::CheckParticleList: " << missingName
             << " do not exist " << G4endl;
      G4cout << " These particle are necessary for basic EM processes"
             << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::CheckParticleList", "Run0101",
                FatalException, "Missing EM basic particle");
  }

  // RULE 2
  //  proton should exist
  //   if any other charged baryon  exist
  if(!isProton && isAnyChargedBaryon)
  {
    G4String missingName = "proton ";

#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::CheckParticleList: " << missingName
             << " does not exist " << G4endl;
      G4cout << " Proton is necessary for EM baryon processes" << G4endl;
    }
#endif
    missingName += " should be created ";
    G4Exception("G4PhysicsListHelper::CheckParticleList", "Run0102",
                FatalException, "Missing Proton");
  }

  // RULE 3
  //  GenericIonn should exist
  //   if any other ion  exist
  if(!isGenericIon && isAnyIon)
  {
    G4String missingName = "GenericIon ";

#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::CheckParticleList: " << missingName
             << " does not exist " << G4endl;
      G4cout << " GenericIon should be created if any ion is necessary"
             << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::CheckParticleList", "Run0103",
                FatalException, "Missing GenericIon");
  }
}

// --------------------------------------------------------------------
void G4PhysicsListHelper::AddTransportation()
{
  G4int verboseLevelTransport = 0;

#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4PhysicsListHelper::AddTransportation()  " << G4endl;
  }
#endif

  G4int nParaWorld =
    G4RunManagerKernel::GetRunManagerKernel()->GetNumberOfParallelWorld();

  if(nParaWorld > 0 || useCoupledTransportation ||
     G4ScoringManager::GetScoringManagerIfExist())
  {
    auto coupledTransport = new G4CoupledTransportation(verboseLevelTransport);
    if(theLooperThresholds == 0)
      coupledTransport->SetLowLooperThresholds();
    if(theLooperThresholds == 2)
      coupledTransport->SetHighLooperThresholds();
    theTransportationProcess = coupledTransport;

    if(verboseLevel > 0)
    {
      G4cout << "--- G4CoupledTransportation is used " << G4endl;
    }
  }
  else
  {
    auto simpleTransport = new G4Transportation(verboseLevelTransport);
    if(theLooperThresholds == 0)
      simpleTransport->SetLowLooperThresholds();
    if(theLooperThresholds == 2)
      simpleTransport->SetHighLooperThresholds();
    theTransportationProcess = simpleTransport;
  }

  // loop over all particles in G4ParticleTable
  aParticleIterator->reset();
  while((*aParticleIterator)())
  {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    // Add transportation process for all particles
    if(pmanager == 0)
    {
      // Error !! no process manager
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4PhysicsListHelper::AddTransportation  "
               << " : No Process Manager for " << particle->GetParticleName()
               << G4endl;
      }
#endif
      G4Exception("G4PhysicsListHelper::AddTransportation", "Run0104",
                  FatalException, "No process manager");
      continue;
    }
    // Molecule use different type transportation
    if(particle->GetParticleType() == "Molecule")
      continue;

    // add transportation with ordering = ( -1, "first", "first" )
    pmanager->AddProcess(theTransportationProcess);
    pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
    pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
  }
}

// --------------------------------------------------------------------
void G4PhysicsListHelper::ReadOrdingParameterTable()
{
  G4bool readInFile = false;
  std::ifstream fIn;

  if(std::getenv("G4ORDPARAMTABLE"))
  {
    ordParamFileName = std::getenv("G4ORDPARAMTABLE");
#ifdef G4VERBOSE
    if(verboseLevel > 1)
    {
      G4cout << "G4PhysicsListHelper::ReadOrdingParameterTable  :"
             << ordParamFileName << " is assigned to Ordering Parameter Table "
             << G4endl;
    }
#endif
    // open input file //
    fIn.open(ordParamFileName, std::ios::in);
    // check if the file has been opened successfully
    if(!fIn)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4PhysicsListHelper::ReadOrdingParameterTable  "
               << " Can not open file " << ordParamFileName << G4endl;
      }
#endif
      G4Exception("G4PhysicsListHelper::ReadOrdingParameterTable", "Run0105",
                  JustWarning, "Fail to open ordering parameter table ");
    }
    else
    {
      readInFile = true;
    }
  }

  // create OrdParamTable
  if(theTable != nullptr)
  {
    theTable->clear();
    delete theTable;
    theTable    = nullptr;
    sizeOfTable = 0;
  }
  theTable    = new G4OrdParamTable();
  sizeOfTable = 0;

  if(readInFile)
  {
    // read in the file and fill the table
    while(!fIn.eof())
    {
      G4PhysicsListOrderingParameter tmp;
      G4int flag;
      fIn >> tmp.processTypeName >> tmp.processType >> tmp.processSubType >>
        tmp.ordering[0] >> tmp.ordering[1] >> tmp.ordering[2] >> flag;
      tmp.isDuplicable = (flag != 0);
      theTable->push_back(tmp);
      ++sizeOfTable;
    }
    fIn.close();
  }
  else
  {
    ReadInDefaultOrderingParameter();
  }

  if(sizeOfTable == 0)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::ReadOrdingParameterTable "
             << " Empty file " << ordParamFileName << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::ReadOrdingParameterTable", "Run0106",
                JustWarning, "The ordering parameter table is empty ");
    delete theTable;
    theTable = nullptr;
  }
  return;
}

// --------------------------------------------------------------------
void G4PhysicsListHelper::DumpOrdingParameterTable(G4int subType) const
{
  if(theTable == nullptr)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::DumpOrdingParameterTable   "
             << " No ordering parameter table  : " << ordParamFileName
             << G4endl;
    }
#endif
    return;
  }
  G4cout << "G4PhysicsListHelper::DumpOrdingParameterTable  : "
         << ordParamFileName << G4endl;
  G4cout << "          TypeName  "
         << "    ProcessType"
         << "        SubType"
         << "         AtRest"
         << "      AlongStep"
         << "        PostStep"
         << "     Duplicable" << G4endl;
  for(G4int i = 0; i < sizeOfTable; ++i)
  {
    G4PhysicsListOrderingParameter* tmp = &(theTable->at(i));
    if((subType >= 0) && (subType != tmp->processSubType))
      continue;
    G4cout << std::setw(18) << tmp->processTypeName << std::setw(15)
           << tmp->processType << std::setw(15) << tmp->processSubType
           << std::setw(15) << tmp->ordering[0] << std::setw(15)
           << tmp->ordering[1] << std::setw(15) << tmp->ordering[2];
    if(tmp->isDuplicable)
    {
      G4cout << "  true";
    }
    else
    {
      G4cout << "  false";
    }
    G4cout << G4endl;
  }
}

// --------------------------------------------------------------------
G4PhysicsListOrderingParameter
G4PhysicsListHelper::GetOrdingParameter(G4int subType) const
{
  G4PhysicsListOrderingParameter value;

  if(theTable == nullptr)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::GetOrderingParameter : "
             << " No ordering parameter table  : " << ordParamFileName
             << G4endl;
    }
#endif
    return value;
  }

  for(G4int i = 0; i < sizeOfTable; ++i)
  {
    G4PhysicsListOrderingParameter* tmp = &(theTable->at(i));
    if(subType == tmp->processSubType)
    {
      value.processTypeName = tmp->processTypeName;
      value.processType     = tmp->processType;
      value.processSubType  = tmp->processSubType;
      value.ordering[0]     = tmp->ordering[0];
      value.ordering[1]     = tmp->ordering[1];
      value.ordering[2]     = tmp->ordering[2];
      value.isDuplicable    = tmp->isDuplicable;
    }
  }
  return value;
}

// --------------------------------------------------------------------
G4bool G4PhysicsListHelper::RegisterProcess(G4VProcess* process,
                                            G4ParticleDefinition* particle)
{
  if(theTable == nullptr)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::RegisterProcess :"
             << " No ordering parameter table  : " << ordParamFileName
             << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterProcess", "Run0107",
                FatalException, "No Ordering Parameter Table");
    return false;
  }

  const G4String pName = process->GetProcessName();
  const G4int pType    = process->GetProcessType();
  const G4int pSubType = process->GetProcessSubType();

#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4PhysicsListHelper::RegisterProcess :" << pName
           << " Process Type = " << pType << " SubType = " << pSubType << " to "
           << particle->GetParticleName() << G4endl;
  }
#endif

  // Check Process Type/SubType
  if((pType < 1) || (pSubType < 1))
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::RegisterProcess :" << pName << " for "
             << particle->GetParticleName()
             << " has illegal Process Type = " << pType
             << " SubType = " << pSubType << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterProcess", "Run0108",
                FatalException, "No Matching process Type/SubType");
    return false;
  }

  G4bool isFound = false;
  G4int ord[3];
  G4bool duplicable = false;
  for(G4int i = 0; i < sizeOfTable; ++i)
  {
    G4PhysicsListOrderingParameter* tmp = &(theTable->at(i));
    if((tmp->processType == pType) && (tmp->processSubType == pSubType))
    {
      ord[0]     = tmp->ordering[0];
      ord[1]     = tmp->ordering[1];
      ord[2]     = tmp->ordering[2];
      duplicable = tmp->isDuplicable;
      isFound    = true;
      break;
    }
  }
  if(!isFound)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::RegisterProcess :" << pName << " for "
             << particle->GetParticleName() << " with  type/subtype =" << pType
             << "/" << pSubType
             << "  is not registered in OrdingParameterTable  " << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterProcess", "Run0109",
                FatalException, "No Matching process Type/SubType");
    return false;
  }

  // Check Process Manager
  G4ProcessManager* pManager = particle->GetProcessManager();
  if(pManager == nullptr)
  {
    // Error !! no process manager
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4PhysicsListHelper::RegisterProcess "
             << " : No Process Manager for " << particle->GetParticleName()
             << G4endl;
    }
#endif
    G4Exception("G4PhysicsListHelper::RegisterProcess   ", "Riun0110",
                FatalException, "No process manager");
    return false;
  }

  // Check Duplication
  if(!duplicable)
  {
    G4bool duplicated      = false;
    G4ProcessVector* pList = pManager->GetProcessList();
    for(G4int idx = 0; idx < (G4int)pList->size(); ++idx)
    {
      const G4VProcess* p = (*pList)[idx];
      if((p->GetProcessType() == pType) && (p->GetProcessSubType() == pSubType))
      {
        duplicated = true;
#ifdef G4VERBOSE
        if(verboseLevel > 0)
        {
          G4cout << "G4PhysicsListHelper::RegisterProcess :" << pName << " for "
                 << particle->GetParticleName()
                 << " with  type/subtype =" << pType << "/" << pSubType
                 << "  is has same subType as " << p->GetProcessName()
                 << " for " << particle->GetParticleName() << G4endl;
          G4cout << "It will not be added !!" << G4endl;
        }
#endif
        G4Exception("G4PhysicsListHelper::RegisterProcess", "Run0111",
                    JustWarning, "Duplication of processes");
      }
    }
    if(duplicated)
      return false;
  }

  // Add Process
  G4int code = pManager->AddProcess(process);
  if(code < 0)
    return false;

  // Set Ordering Parameter
  for(G4int idx = 0; idx < 3; ++idx)
  {
    G4ProcessVectorDoItIndex idxOrd =
      static_cast<G4ProcessVectorDoItIndex>(idx);
    if(ord[idx] < 0)
    {
      // Do Nothing because NO DOIT
    }
    else if(ord[idx] == 0)
    {
      pManager->SetProcessOrderingToFirst(process, idxOrd);
    }
    else if(ord[idx] < 9999)
    {
      pManager->SetProcessOrdering(process, idxOrd, ord[idx]);
    }
    else
    {
      pManager->SetProcessOrderingToLast(process, idxOrd);
    }
  }
#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "G4PhysicsListHelper::RegisterProcess :" << pName << " for "
           << particle->GetParticleName() << " with  type/subtype =" << pType
           << "/" << pSubType
           << " is successfully registered with ordering parameters " << ord[0]
           << ":" << ord[1] << ":" << ord[2] << G4endl;
  }
#endif
  return true;
}

// --------------------------------------------------------------------
void G4PhysicsListHelper::ReadInDefaultOrderingParameter()
{
  G4PhysicsListOrderingParameter tmp;

  // NOTE: please use enum values, rather than numerical values, 
  //       for both the processType and processSubType below.

  tmp.processTypeName = "Transportation";
  tmp.processType     = fTransportation;
  tmp.processSubType  = TRANSPORTATION;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 0;
  tmp.ordering[2]     = 0;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "CoupleTrans";
  tmp.processType     = fTransportation;
  tmp.processSubType  = COUPLED_TRANSPORTATION;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 0;
  tmp.ordering[2]     = 0;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "CoulombScat";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fCoulombScattering;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Ionisation";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fIonisation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 2;
  tmp.ordering[2]     = 2;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Brems";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fBremsstrahlung;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 3;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "PairProdCharged";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fPairProdByCharged;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 4;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Annih";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fAnnihilation;
  tmp.ordering[0]     = 5;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 5;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "AnnihToMuMu";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fAnnihilationToMuMu;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 6;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "AnnihToTauTau";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fAnnihilationToTauTau;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 7;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "AnnihToHad";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fAnnihilationToHadrons;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 8;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "NuclearStopping";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fNuclearStopping;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 9;
  tmp.ordering[2]     = -1;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "ElectronGeneral";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fElectronGeneralProcess;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 1;
  tmp.ordering[2]     = 1;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Msc";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fMultipleScattering;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 1;
  tmp.ordering[2]     = -1;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Rayleigh";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fRayleigh;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "PhotoElectric";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fPhotoElectricEffect;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Compton";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fComptonScattering;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Conv";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fGammaConversion;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "ConvToMuMu";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fGammaConversionToMuMu;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "GammaGeneral";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fGammaGeneralProcess;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "PositronGeneral";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fPositronGeneralProcess;
  tmp.ordering[0]     = 1;
  tmp.ordering[1]     = 1;
  tmp.ordering[2]     = 1;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "MuPairByMuon";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fMuonPairProdByCharged;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 10;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Cerenkov";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fCerenkov;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Scintillation";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fScintillation;
  tmp.ordering[0]     = 9999;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 9999;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "SynchRad";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fSynchrotronRadiation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "TransRad";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fTransitionRadiation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "SurfaceRefl";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fSurfaceReflection;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "OpAbsorb";
  tmp.processType     = fOptical;
  tmp.processSubType  = fOpAbsorption;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "OpBoundary";
  tmp.processType     = fOptical;
  tmp.processSubType  = fOpBoundary;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "OpRayleigh";
  tmp.processType     = fOptical;
  tmp.processSubType  = fOpRayleigh;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "OpWLS";
  tmp.processType     = fOptical;
  tmp.processSubType  = fOpWLS;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "OpMieHG";
  tmp.processType     = fOptical;
  tmp.processSubType  = fOpMieHG;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "OpWLS2";
  tmp.processType     = fOptical;
  tmp.processSubType  = fOpWLS2;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAElastic";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyElastic;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAExcit";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyExcitation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAIonisation";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyIonisation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAVibExcit";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyVibrationalExcitation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAAttachment";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyAttachment;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAChargeDec";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyChargeDecrease;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAChargeInc";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyChargeIncrease;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAElecSolv";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyElectronSolvation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAMolecDecay";
  tmp.processType     = fDecay;
  tmp.processSubType  = fLowEnergyMolecularDecay;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = -1;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "ITTransport";
  tmp.processType     = fTransportation;
  tmp.processSubType  = fLowEnergyTransportation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 0;
  tmp.ordering[2]     = 0;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNABrownTrans";
  tmp.processType     = fTransportation;
  tmp.processSubType  = fLowEnergyBrownianTransportation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = 0;
  tmp.ordering[2]     = 0;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNADoubleIoni";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyDoubleIonisation;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNADoubleCap";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyDoubleCap;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAIoniTransfer";
  tmp.processType     = fElectromagnetic;
  tmp.processSubType  = fLowEnergyIoniTransfer;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DNAStaticMol";
  tmp.processType     = fUserDefined;
  tmp.processSubType  = fLowEnergyStaticMol;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     =  1000;
  tmp.isDuplicable =  false;
  theTable->push_back(tmp);
  sizeOfTable +=1;

  tmp.processTypeName = "DNAScavenger";
  tmp.processType     = fUserDefined;
  tmp.processSubType  = fLowEnergyScavenger;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "HadElastic";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fHadronElastic;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "NeutronGeneral";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fNeutronGeneral;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "HadInelastic";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fHadronInelastic;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "HadCapture";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fCapture;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "MuAtomCapture";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fMuAtomicCapture;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "HadFission";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fFission;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "HadAtRest";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fHadronAtRest;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = -1;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "HadCEX";
  tmp.processType     = fHadronic;
  tmp.processSubType  = fChargeExchange;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "Decay";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DecayWSpin";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY_WithSpin;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DecayPiSpin";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY_PionMakeSpin;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DecayRadio";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY_Radioactive;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DecayUnKnown";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY_Unknown;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DecayMuAtom";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY_MuAtom;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "DecayExt";
  tmp.processType     = fDecay;
  tmp.processSubType  = DECAY_External;
  tmp.ordering[0]     = 1000;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "StepLimiter";
  tmp.processType     = fGeneral;
  tmp.processSubType  = STEP_LIMITER;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "UsrSepcCuts";
  tmp.processType     = fGeneral;
  tmp.processSubType  = USER_SPECIAL_CUTS;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "NeutronKiller";
  tmp.processType     = fGeneral;
  tmp.processSubType  = NEUTRON_KILLER;
  tmp.ordering[0]     = -1;
  tmp.ordering[1]     = -1;
  tmp.ordering[2]     = 1000;
  tmp.isDuplicable    = false;
  theTable->push_back(tmp);
  sizeOfTable += 1;

  tmp.processTypeName = "ParallelWorld";
  tmp.processType     = fParallel;
  tmp.processSubType  = PARALLEL_WORLD_PROCESS;
  tmp.ordering[0]     = 9900;
  tmp.ordering[1]     = 1;
  tmp.ordering[2]     = 9900;
  tmp.isDuplicable    = true;
  theTable->push_back(tmp);
  sizeOfTable += 1;
}
