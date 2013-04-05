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

namespace  {
    //This is used to lock on shared resource
    G4Mutex xsFactoryRegisterMutex = G4MUTEX_INITIALIZER;
}

G4CrossSectionFactoryRegistry* G4CrossSectionFactoryRegistry::instance = 0;
G4CrossSectionFactoryRegistry* G4CrossSectionFactoryRegistry::Instance()
{
    G4AutoLock l(&xsFactoryRegisterMutex);
    if ( !instance ) instance = new G4CrossSectionFactoryRegistry();
    return instance;
}

G4CrossSectionFactoryRegistry::G4CrossSectionFactoryRegistry()
{
}

G4CrossSectionFactoryRegistry::G4CrossSectionFactoryRegistry(const G4CrossSectionFactoryRegistry&)
{
    G4Exception("G4CrossSectionFactoryRegistry::G4CrossSectionFactoryRegistry",
                "CrossSection004",FatalException,"Use of copy constructor not allowed");
}
G4CrossSectionFactoryRegistry& G4CrossSectionFactoryRegistry::operator=(const G4CrossSectionFactoryRegistry&)
{
    G4Exception("G4CrossSectionFactoryRegistry::G4CrossSectionFactoryRegistry",
                "CrossSection004",FatalException,"Use of assignement operator not allowed");
    return *this;
}

void G4CrossSectionFactoryRegistry::Register( const G4String& name, G4VBaseXSFactory* factory )
{
    G4AutoLock l(&xsFactoryRegisterMutex);
    if ( factories.find(name) != factories.end() )
    {
        G4ExceptionDescription msg;
        msg <<"Cross section factory with name: "<<name
            <<" already existing, old factory has been replaced";
        G4Exception("G4CrossSectionFactoryRegistry::Register(...)",
                    "CrossSection002",JustWarning,msg);
    }
    factories[name] = factory;
}

G4VBaseXSFactory* G4CrossSectionFactoryRegistry::GetFactory( const G4String& name, G4bool abortIfNotFound ) const
{
    G4AutoLock l(&xsFactoryRegisterMutex);
    std::map<G4String,G4VBaseXSFactory*>::const_iterator it = factories.find(name);
    if ( it != factories.end() ) return it->second;
    else
    {
        if ( abortIfNotFound )
        {
            G4ExceptionDescription msg;
            msg <<"Cross section factory with name: "<<name
                <<" not found.";
            G4Exception("G4CrossSectionFactoryRegistry::Register(...)",
                        "CrossSection003",FatalException,msg);
        }
    }
    return static_cast<G4VBaseXSFactory*>(0);
}


std::ostream&  operator<<(std::ostream& msg, const G4CrossSectionFactoryRegistry& rhs) {
    msg<<"Factory Registry "<<&rhs<<" has factories: [";
    for ( std::map<G4String,G4VBaseXSFactory*>::const_iterator it =rhs.factories.begin() ;
         it != rhs.factories.end() ; ++it )
    {
        msg<<(*it).first<<":"<<(*it).second<<",";
    }
    msg<<"]";
    return msg;
}
