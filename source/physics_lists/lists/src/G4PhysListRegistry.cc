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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:    G4PhysListRegistry
//
// Author  R. Hatcher  2014-10-15
//
// Modifications:  based on G4PhysicsConstructorRegistry
//
//#define G4VERBOSE 1

#include "G4ios.hh"
#include <iomanip>
#include <algorithm>

#include "G4PhysListRegistry.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PhysListStamper.hh"

#include "G4PhysicsConstructorRegistry.hh"
// include G4_REFERENCE_PHYSCONSTR_FACTORY() macros here to pull them in
// for static builds.   Would do this in  G4PhysicsConstructorRegistry
// but that causes a circular dependency when the phys_ctor libs
// are broken in a "granular" build.
#define REGREF PhysListReg // REGREF is used to make instances unique
#include "G4RegisterPhysicsConstructors.icc"

// include this with this compilation unit so static builds pull them in
#include "G4RegisterPhysLists.icc"

G4ThreadLocal G4PhysListRegistry* G4PhysListRegistry::theInstance = 0;

const int doReplace  = 0x01;  // _ (ReplacePhysics) vs. + (RegisterPhysics)
const int isCtorName = 0x02;  // true if actual physCtor name vs. ext name

G4PhysListRegistry* G4PhysListRegistry::Instance()
{
  if ( 0 == theInstance) {
    static G4ThreadLocal G4PhysListRegistry *manager_G4MT_TLS_ = 0;
    if (!manager_G4MT_TLS_) manager_G4MT_TLS_ = new G4PhysListRegistry;
    G4PhysListRegistry &manager = *manager_G4MT_TLS_;
    theInstance = &manager;
  }

  // common EM overrides
  theInstance->AddPhysicsExtension("EM0","G4EmStandardPhysics");
  theInstance->AddPhysicsExtension("EMV","G4EmStandardPhysics_option1");
  theInstance->AddPhysicsExtension("EMX","G4EmStandardPhysics_option2");
  theInstance->AddPhysicsExtension("EMY","G4EmStandardPhysics_option3");
  theInstance->AddPhysicsExtension("EMZ","G4EmStandardPhysics_option4");
  theInstance->AddPhysicsExtension("LIV","G4EmLivermorePhysics");
  theInstance->AddPhysicsExtension("PEN","G4EmPenelopePhysics");
  // the GS EM extension originally required double underscores
  // support either one or two as __GS is confusing to users
  // same for __SS
  theInstance->AddPhysicsExtension("GS" ,"G4EmStandardPhysicsGS");
  theInstance->AddPhysicsExtension("_GS","G4EmStandardPhysicsGS");
  theInstance->AddPhysicsExtension("SS" ,"G4EmStandardPhysicsSS");
  theInstance->AddPhysicsExtension("_SS","G4EmStandardPhysicsSS");

  theInstance->AddPhysicsExtension("EM0","G4EmStandardPhysics");
  theInstance->AddPhysicsExtension("WVI","G4EmStandardPhysicsWVI");
  theInstance->AddPhysicsExtension("LE" ,"G4EmLowEPPhysics");
  theInstance->AddPhysicsExtension("_LE","G4EmLowEPPhysics");

  return theInstance;
}

G4PhysListRegistry::G4PhysListRegistry()
  : verbose(1)
  , unknownFatal(0)
  , systemDefault("FTFP_BERT")
{
  SetUserDefaultPhysList();
}

G4PhysListRegistry::~G4PhysListRegistry()
{
}

void G4PhysListRegistry::SetUserDefaultPhysList(const G4String& name)
{
  if ( name == "" ) userDefault = systemDefault;
  else              userDefault = name;
}

void G4PhysListRegistry::AddFactory(G4String name, G4VBasePhysListStamper* factory)
{
  factories[name] = factory;
}

void G4PhysListRegistry::AddPhysicsExtension(G4String name, G4String procname)
{
  // a mapping from short extension names to actual physics process constructors
  physicsExtensions[name] = procname;
}

G4VModularPhysicsList*
G4PhysListRegistry::GetModularPhysicsList(const G4String& name)
{
  //
  //
  G4String plBase = "";
  std::vector<G4String> physExt;
  std::vector<G4int>    physReplace;
  G4bool allKnown =
    DeconstructPhysListName(name,plBase,physExt,physReplace,verbose);

  size_t npc = physExt.size();
  if ( verbose > 0 ) {
    G4cout << "G4PhysListRegistry::GetModularPhysicsList <"
           << name << ">"
           << ", as \"" << plBase << "\" with extensions \"";
    for ( size_t ipc = 0; ipc < npc; ++ipc )
      G4cout << ((physReplace[ipc]&doReplace)?"_":"+") << physExt[ipc];
    G4cout << "\"" << G4endl;
  }

  if ( ! allKnown ) {
    // couldn't match what the user wanted ...
    G4cout << "### G4PhysListRegistry WARNING: " << name
           << " is not known" << G4endl << G4endl;
    if ( ! unknownFatal ) return 0;

    G4ExceptionDescription ED;
    ED << "The factory for the physicslist ["<< name << "] does not exist!"
       << G4endl;
    if ( plBase == "" ) {
      ED << "Could determine no sensible base physics list" << G4endl;
    } else {
      ED << "One or more of the extensions does not exist [ ";
      for ( size_t ipc = 0; ipc < physExt.size(); ++ipc ) {
        ED << physExt[ipc] << " ";
      }
      ED << "]" << G4endl;
    }
    G4Exception("G4PhysListRegistry::GetModularPhysicsList",
                "PhysicsList002", FatalException, ED);
    return 0;
  }

  // if we want this method "const" then the next line becomes more complex
  // because there is no const version of [] (which adds an entry if the
  // key doesn't exist)
  G4VModularPhysicsList* pl = factories[plBase]->Instantiate(verbose);
  G4PhysicsConstructorRegistry* pcRegistry =
    G4PhysicsConstructorRegistry::Instance();
  G4int ver = pl->GetVerboseLevel();
  pl->SetVerboseLevel(0);
  for ( size_t ipc = 0; ipc < npc; ++ipc ) {
    // got back a list of short names, need to use the map to get the
    // full physics constructor name
    G4String extName = physExt[ipc];
    G4String pcname  =
      ((physReplace[ipc]&isCtorName)) ? extName : physicsExtensions[extName];
    // this doesn't have a verbose option ... it should
    // but G4PhysicsConstructorFactory doesn't support it
    G4VPhysicsConstructor* pctor = pcRegistry->GetPhysicsConstructor(pcname);
    G4String reporreg = "";
    if (( physReplace[ipc] & doReplace)) {
      pl->ReplacePhysics(pctor);
      reporreg = "ReplacePhysics ";
    } else {
      pl->RegisterPhysics(pctor);
      reporreg = "RegisterPhysics";
    }
    if ( verbose > 0 ) G4cout << "<<< " << reporreg << " with " << pcname
                              << " \"" << extName << "\"" << G4endl;
  }
  pl->SetVerboseLevel(ver);
  G4cout << "<<< Reference Physics List " << name << " is built" << G4endl;
  G4cout << G4endl; // old factory has this

  return pl;
}

G4VModularPhysicsList*
G4PhysListRegistry::GetModularPhysicsListFromEnv()
{
  //
  // instantiate PhysList by environment variable "PHYSLIST"
  // if not set use default
  G4String name = "";
  char* path = std::getenv("PHYSLIST");
  if (path) {
    name = G4String(path);
  } else {
    name = userDefault;
    G4cout << "### G4PhysListRegistry WARNING: "
           << " environment variable PHYSLIST is not defined"
           << G4endl
           << "    Default Physics Lists " << name
           << " is instantiated"
           << G4endl;
  }
  return GetModularPhysicsList(name);
}

G4bool G4PhysListRegistry::IsReferencePhysList(G4String name) const
{
  G4String plBase = "";
  std::vector<G4String> physExt;
  std::vector<G4int>    physReplace;
  G4bool allKnown = DeconstructPhysListName(name,plBase,physExt,physReplace,1);
  return allKnown;
}

G4bool G4PhysListRegistry::DeconstructPhysListName(const G4String& name,
                                                   G4String& plBase,
                                                   std::vector<G4String>& physExt,
                                                   std::vector<G4int>& replace,
                                                   G4int verb) const
{
  // Take apart a name given to us by the user
  // this name might be a base PhysList + unknown number of extensions
  // Extensions preceeded with a "_" should use
  //    ReplacePhysics()
  // those proceeded with a "+" should use
  //    RegisterPhysics()
  // the former is in line with previous behaviour, while the second allows
  // additional flexibility
  plBase = "";
  physExt.clear();
  replace.clear();
  bool allKnown = false;

  G4String workingName = name;

  const std::vector<G4String>& availBases  = AvailablePhysLists();
  const std::vector<G4String>& availExtras = AvailablePhysicsExtensions();

  G4PhysicsConstructorRegistry* physConstRegistry =
    G4PhysicsConstructorRegistry::Instance();

  const std::vector<G4String>& availPhysCtors  =
    physConstRegistry->AvailablePhysicsConstructors();

  // find the longest base list that is contained in the user supplied name
  // and starts at the beginning
  G4String bestBase = "";
  allKnown = FindLongestMatch(workingName,"base",availBases,plBase);
  if ( verb > 2 ) {
    G4cout << "  " << name << ", base known=" << ((allKnown)?"true":"false")
           << " chosen plBase \"" << plBase << "\"" << G4endl;
  }
  if ( ! allKnown ) {
    // didn't find any matching base physics list
    // no point of going on to the extensions
    return allKnown;
  }
  // remove base name for working name
  workingName.erase(0,plBase.size());

  // now start trying to match up extensions and/or physCtors
  // each should be preceeded by at "_" (replace) or "+" (register)
  // but don't freak if it isn't, just assume "_"
  while ( ! workingName.empty() ) {
    char c = workingName.data()[0];  // leading character
    if ( '_' == c || '+' == c ) workingName.erase(0,1);  // and remove it
    G4int    replaceExtra = (( c != '+' ) ? doReplace : 0 );
    G4String extraName = "";
    G4bool   extraKnown = false;

    extraKnown = FindLongestMatch(workingName,"extNames",availExtras,extraName);
    if ( extraKnown ) {
      // physics mapping name is known, but is it actually linked to physics?
      //const issue// G4String pcname = physicsExtensions[extraName];
      std::map<G4String,G4String>::const_iterator itr =
        physicsExtensions.find(extraName);
      G4String pcname = "";
      if ( itr != physicsExtensions.end() ) pcname = itr->second;
      bool realknown = physConstRegistry->IsKnownPhysicsConstructor(pcname);
      if ( ! realknown ) allKnown = false;
#ifdef G4VERBOSE
      if ( verb > 2 ) {
        G4cout << "  extraName \"" << extraName << "\" maps to physics ctor \""
               << pcname << "\" which is itself realknown " << realknown
               << G4endl;
      }
#endif
    } else {
      // perhaps it's an explicit physCtor name
      extraKnown =
        FindLongestMatch(workingName,"physCtors",availPhysCtors,extraName);
      if ( extraKnown ) replaceExtra |= isCtorName; // flag it
    }
#ifdef G4VERBOSE
    if ( verb > 2 ) {
      G4cout << "  physextra " << name << " [" << workingName << "]"
             <<", extra known " << extraKnown
             << " chosen extra \"" << extraName << "\""
             << " replace " << replaceExtra << G4endl;
    }
#endif
    if ( extraKnown ) {
      physExt.push_back(extraName);
      replace.push_back(replaceExtra);
      // and remove it so we can look for the next bit
      workingName.erase(0,extraName.size());

    } else {
#ifdef G4VERBOSE
      if ( verb > 2 ) {
        G4cout << "  workingName \"" << workingName << "\""
               << " couldn't be found in the extensions list"
               << G4endl;
      }
#endif
      allKnown = false;
      // found a pattern that we can't map
      return allKnown;
    }
  } // workingName not empty

  return allKnown;
}

G4bool G4PhysListRegistry::FindLongestMatch(const G4String& workingName,
                                            const G4String& searchName,
                                            const std::vector<G4String>& validNames,
                                            G4String& bestMatch,
                                            G4int verb) const
{
  bestMatch = "";
  bool found = false;

  size_t n = validNames.size();
  for (size_t i=0;  i<n; ++i) {
    const G4String& testName = validNames[i];
    size_t ipos = workingName.find(testName);
    if ( ipos == 0 ) {
      if ( testName.size() > bestMatch.size() )  {
        bestMatch = testName;
        found     = true;
        if ( verb > 3 ) {
          G4cout << "  " << searchName << " current best guess: "
                 << testName << G4endl;
        }
      } else {
        if ( verb > 3 ) {
          G4cout << "  " << searchName << " match but shorter: "
                 << testName << G4endl;
        }
      }
    } else {
      if ( verb > 3 ) {
        G4cout << "  " << searchName << " reject: " << testName << G4endl;
      }
    }
  }
  return found;
}


const std::vector<G4String>& G4PhysListRegistry::AvailablePhysLists() const
{
  availBasePhysLists.clear();
  std::map<G4String,G4VBasePhysListStamper*>::const_iterator itr;
  for ( itr = factories.begin(); itr != factories.end(); ++itr ) {
    availBasePhysLists.push_back(itr->first);
  }

  return availBasePhysLists;
}

const std::vector<G4String>& G4PhysListRegistry::AvailablePhysicsExtensions() const
{
  availExtensions.clear();
  std::map<G4String,G4String>::const_iterator itr;
  for ( itr = physicsExtensions.begin(); itr != physicsExtensions.end(); ++itr ) {
    availExtensions.push_back(itr->first);
  }

  return availExtensions;
}

const std::vector<G4String>& G4PhysListRegistry::AvailablePhysListsEM() const
{
  // in principle this method could weed out all the extensions that aren't
  // EM replacements ... but for now just use it as a synonym for
  // AvailablePhysicsExtensions()
  return AvailablePhysicsExtensions();
}

void G4PhysListRegistry::PrintAvailablePhysLists() const
{
  std::vector<G4String> avail = AvailablePhysLists();
  G4cout << "Base G4VModularPhysicsLists in G4PhysListRegistry are:"
         << G4endl;
  if ( avail.empty() ) G4cout << "... no registered lists" << G4endl;
  else {
    size_t n = avail.size();
    for (size_t i=0; i<n; ++i ) {
      G4cout << " [" << std::setw(3) << i << "] "
             << " \"" << avail[i] << "\"" << G4endl;
    }
  }

  G4PhysicsConstructorRegistry* physConstRegistry = G4PhysicsConstructorRegistry::Instance();

  std::map<G4String,G4String>::const_iterator itr;
  G4cout << "Replacement mappings in G4PhysListRegistry are:"
         << G4endl;
  for ( itr = physicsExtensions.begin(); itr != physicsExtensions.end(); ++itr ) {
    bool known = physConstRegistry->IsKnownPhysicsConstructor(itr->second);

    G4cout << "    " << std::setw(10) << itr->first << " => "
           << std::setw(30) << itr->second << " "
           << ( (known)?"":"[unregistered physics]")
           << G4endl;
  }
  G4cout << "Use these mapping to extend physics list; append with _EXT or +EXT" << G4endl
         << "   to use ReplacePhysics() (\"_\") or RegisterPhysics() (\"+\")."
         << G4endl;
}
