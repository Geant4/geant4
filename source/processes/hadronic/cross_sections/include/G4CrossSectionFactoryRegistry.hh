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
// This class implements a regsitry for cross-section factories
// It can be used to access factories to instantiate cross-section
// datasets on demand. This is intended to be shared among threads
// see G4CrossSectionFactory class: while G4CrossSectionDataSetRegistry
// is per-thread (G4ThreadLocal), this class is shared, so only one instance
// of each factory exists for each application. This is needed because factories
// rely on global variables to register themselves in this registry.
// Note: use of this class is thread-safe.
//
// History:
// 1-Apr-2013: A. Dotti, first implementation
#ifndef G4CROSSSECTIONFACTORYREGISTRY_HH
#define G4CROSSSECTIONFACTORYREGISTRY_HH

#include <ostream>
#include <G4String.hh>
#include <map>

class G4VBaseXSFactory;

class G4CrossSectionFactoryRegistry
{
    friend std::ostream& operator<<(std::ostream&, const G4CrossSectionFactoryRegistry&);
private:
    std::map<G4String, G4VBaseXSFactory*> factories;
    static G4CrossSectionFactoryRegistry* instance; //Note this is shared among threads
    G4CrossSectionFactoryRegistry();
    G4CrossSectionFactoryRegistry(const G4CrossSectionFactoryRegistry& );
    G4CrossSectionFactoryRegistry& operator=(const G4CrossSectionFactoryRegistry&);
    //Disable copy-ctr and assignement operator
public:
    static G4CrossSectionFactoryRegistry* Instance();
    G4VBaseXSFactory* GetFactory( const G4String& name , G4bool abortIfNotFound = true) const;
    //Search a cross-section factory by name, by default rise an exception if factory is not found 
    void Register( const G4String& name , G4VBaseXSFactory* factory );
};

std::ostream&  operator<<(std::ostream& msg, const G4CrossSectionFactoryRegistry& rhs);

#endif // G4CROSSSECTIONFACTORYREGISTRY_HH
