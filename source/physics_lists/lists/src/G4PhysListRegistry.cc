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
// $Id: 
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4PhysListRegistry
//
// Author  R. Hatcher  2014-10-15
//
// Modifications:  based on G4PhysicsConstructorRegistry
//

#include "G4ios.hh"
#include <iomanip>
#include <algorithm>

#include "G4PhysListRegistry.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PhysListStamper.hh"
#include "G4PhysicsConstructorRegistry.hh"

G4ThreadLocal G4PhysListRegistry* G4PhysListRegistry::theInstance = 0;

G4PhysListRegistry* G4PhysListRegistry::Instance()
{
  if ( 0 == theInstance) {
    static G4ThreadLocal G4PhysListRegistry *manager_G4MT_TLS_ = 0; 
    if (!manager_G4MT_TLS_) manager_G4MT_TLS_ = new G4PhysListRegistry;
    G4PhysListRegistry &manager = *manager_G4MT_TLS_;
    theInstance = &manager;
  }

  // common EM overrides
  theInstance->AddPhysicsExtension("EMV","G4EmStandardPhysics_option1");
  theInstance->AddPhysicsExtension("EMX","G4EmStandardPhysics_option2");
  theInstance->AddPhysicsExtension("EMY","G4EmStandardPhysics_option3");
  theInstance->AddPhysicsExtension("EMZ","G4EmStandardPhysics_option4");
  theInstance->AddPhysicsExtension("LIV","G4EmLivermorePhysics");
  theInstance->AddPhysicsExtension("PEN","G4EmPenelopePhysics");

  return theInstance;
}

G4PhysListRegistry::G4PhysListRegistry() : verbose(1), unknownFatal(0)
{}

G4PhysListRegistry::~G4PhysListRegistry()
{
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
  G4bool allKnown = DeconstructPhysListName(name,plBase,physExt,physReplace,1);

  if ( ! allKnown ) {
    // couldn't match what the user wanted ...
    if ( ! unknownFatal ) return 0;

    G4ExceptionDescription ED;
    ED << "The factory for the physicslist ["<< name << "] does not exist!" << G4endl;
    if ( plBase == "" ) {
      ED << "Could determine no sensible base physics list" << G4endl;
    } else {
      ED << "One or more of the extensions does not exist [ ";
      for ( size_t ipc = 0; ipc < physExt.size(); ++ipc ) {
        ED << physExt[ipc] << " ";
      }
      ED << "]" << G4endl;
    }
    G4Exception("G4PhysListRegistry::GetModularPhysicsList", "PhysicsList001", FatalException, ED);
    return 0;
  }

  // if we want this method "const" then the next line becomes more complex
  // because there is no const version of [] (which adds an entry if the key doesn't exist)
  G4VModularPhysicsList* pl = factories[plBase]->Instantiate(verbose);
  G4PhysicsConstructorRegistry* pcRegistry = G4PhysicsConstructorRegistry::Instance();
  size_t npc = physExt.size();
  for ( size_t ipc = 0; ipc < npc; ++ipc ) {
    // got back a list of short names, need to use the map to get the full physics constructor name
    G4String extName = physExt[ipc];
    G4String pcname = physicsExtensions[extName];
    G4VPhysicsConstructor* pc = pcRegistry->GetPhysicsConstructor(pcname);
    if ( physReplace[ipc] > 0 ) {
      pl->ReplacePhysics(pc);
    } else {
      pl->RegisterPhysics(pc);
    }
  }
  return pl;
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

  const std::vector<G4String>& availBases = AvailablePhysLists();
  const std::vector<G4String>& availExtras = AvailablePhysicsExtensions();

  // find the longest base list that is contained in the user supplied name
  // and starts at the beginning
  size_t nb = availBases.size();
  for (size_t ib=0;  ib<nb; ++ib) {
    const G4String& testBase = availBases[ib];
    size_t ipos = workingName.find(testBase);
    if ( ipos == 0 ) {
      if ( testBase.size() > plBase.size() )  {
        plBase = testBase;
        allKnown = true;
        if ( verb > 3 ) { G4cout << "  physlist current best guess: " << testBase << G4endl; }
      } else {
        if ( verb > 3 ) { G4cout << "  physlist match but shorter: " << testBase << G4endl; }
      }
    } else {
      if ( verb > 3 ) { G4cout << "  physlist reject: " << testBase << G4endl; }
    }
  }
  if ( verb > 2 ) { 
    G4cout << "  physlist " << name << ", base known " << allKnown 
           << " chosen plBase \"" << plBase << "\"" << G4endl; 
  }
  if ( ! allKnown ) {
    // didn't find any matching base physics list
    // no point of going on to the extensions
    return allKnown;
  }
  // remove base name for working name
  workingName.erase(0,plBase.size());

  // now start trying to match up extensions
  // each should be preceeded by at "_" (replace) or "+" (register)
  // but don't freak if it isn't, just assume "_"
  size_t ne = availExtras.size();
  while ( ! workingName.empty() ) {
    char c = workingName.data()[0];  // leading character
    if ( '_' == c || '+' == c ) workingName.erase(0,1);  // and remove it
    G4int    replaceExtra = ( c != '+' );
    G4String extraName = "";
    G4bool   extraKnown = false;
    for (size_t ie=0;  ie<ne; ++ie) {
      const G4String& testExtra = availExtras[ie];
      size_t ipos = workingName.find(testExtra);
      if ( ipos == 0 ) {
        if ( testExtra.size() > extraName.size() )  {
          extraName = testExtra;
          extraKnown = true;
#ifdef G4VERBOSE
          if ( verb > 3 ) { G4cout << "  physextra current best guess: " << testExtra << G4endl; }
        } else {
          if ( verb > 3 ) { G4cout << "  physextra match but shorter: " << testExtra << G4endl; }
#endif
        }
      } else {
#ifdef G4VERBOSE
        if ( verb > 3 ) { G4cout << "  physextra reject: " << testExtra << G4endl; }
#endif
      }
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
      allKnown = false;
    }
  } // workingName not empty

  return allKnown;
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
  // EM replacements ... but for now just use it as a synonym for AvailablePhysicsExtensions()
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
  ///std::vector<G4String> physConstNames = physConstRegistry->AvailablePhysicsConstructors();

  std::map<G4String,G4String>::const_iterator itr;
  G4cout << "Replacement mappings in G4PhysListRegistry are:"
         << G4endl;
  for ( itr = physicsExtensions.begin(); itr != physicsExtensions.end(); ++itr ) {
    //bool known = ( find(physConstNames.begin(),physConstName.end(),itr->second) != physConstNames.end() );
    bool known = physConstRegistry->IsKnownPhysicsConstructor(itr->second);

    G4cout << "    " << std::setw(10) << itr->first << " => "
           << std::setw(30) << itr->second << " "
           << ( (known)?"":"[unregistered physics]")
           << G4endl;

  }

}

// need something that will trigger compilation
//#include "RegisterPhysLists.icc"

