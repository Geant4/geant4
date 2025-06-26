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
// G4CrossSectionFactoryRegistry class
//
// History:
// 1-Apr-2013: A. Dotti, first implementation
#include "G4CrossSectionFactoryRegistry.hh"
#include "G4CrossSectionFactory.hh"
#include "G4AutoLock.hh"
#include "globals.hh"

namespace
{
  G4Mutex regMutex = G4MUTEX_INITIALIZER;
}

G4CrossSectionFactoryRegistry* G4CrossSectionFactoryRegistry::Instance()
{
  static G4CrossSectionFactoryRegistry reg;
  return &reg;
}

G4CrossSectionFactoryRegistry::G4CrossSectionFactoryRegistry()
{}

void G4CrossSectionFactoryRegistry::Register( const G4String& name, G4VBaseXSFactory* factory )
{
  if ( factories.find(name) == factories.end() ) {
    G4AutoLock l(regMutex);
    if ( factories.find(name) == factories.end() ) {
      factories[name] = factory;
    }
    l.unlock();
  }
}

void G4CrossSectionFactoryRegistry::DeRegister( G4VBaseXSFactory* factory )
{
  if ( nullptr == factory || factories.empty() ) { return; }
  G4AutoLock l(regMutex);
  for ( auto const & f : factories ) {
    if ( factory == f.second ) {
      factories[f.first] = nullptr;
    }
  }
  l.unlock();
}

G4VBaseXSFactory*
G4CrossSectionFactoryRegistry::GetFactory( const G4String& name, G4bool abortIfNotFound ) const
{
  G4VBaseXSFactory* ptr = nullptr;
  auto it = factories.find(name);
  if ( it != factories.end() ) {
    ptr = it->second;
  } else if ( abortIfNotFound ) {
    G4ExceptionDescription msg;
    msg << "Cross section factory with name: " << name << " not found.";
    G4Exception("G4CrossSectionFactoryRegistry::GetFactory(...)",
		"CrossSection003", FatalException, msg);
  }
  return ptr;
}

std::ostream&  operator<<(std::ostream& msg, const G4CrossSectionFactoryRegistry& rhs) {
  msg<<"Factory Registry "<<&rhs<<" has factories: [";
  for ( std::map<G4String,G4VBaseXSFactory*>::const_iterator it =rhs.factories.begin();
	it != rhs.factories.end() ; ++it )
    {
        msg<<(*it).first<<":"<<(*it).second<<",";
    }
    msg<<"]";
    return msg;
}
