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
// G4HCIOcatalog implementation
//
// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------

#include "G4HCIOcatalog.hh"
#include "G4VHCIOentry.hh"

G4ThreadLocal G4HCIOcatalog* G4HCIOcatalog::f_thePointer = 0;

// --------------------------------------------------------------------
G4HCIOcatalog::G4HCIOcatalog()
{
}

// --------------------------------------------------------------------
G4HCIOcatalog* G4HCIOcatalog::GetHCIOcatalog()
{
  if(f_thePointer == nullptr)
    f_thePointer = new G4HCIOcatalog;
  return f_thePointer;
}

// --------------------------------------------------------------------
void G4HCIOcatalog::RegisterEntry(G4VHCIOentry* d)
{
  if(m_verbose > 0)
  {
    G4cout << "registering I/O manager entry \"" << d->GetName() << "\" " << d
           << "." << G4endl;
  }
  if(theCatalog.find(d->GetName()) != theCatalog.cend())
  {
    G4cout << "Redefining I/O Managers list " << d->GetName() << G4endl;
  }
  else
  {
    theCatalog[d->GetName()] = d;
  }
}

// --------------------------------------------------------------------
void G4HCIOcatalog::RegisterHCIOmanager(G4VPHitsCollectionIO* d)
{
  if(m_verbose > 0)
  {
    G4cout << "registering I/O manager \"" << d->SDname() << "\" " << d << "."
           << G4endl;
  }
  if(theStore.find(d->SDname()) != theStore.cend())
  {
    G4cout << "Redefining I/O Manager " << d->SDname() << G4endl;
  }
  else
  {
    theStore[d->SDname()] = d;
  }
}

// --------------------------------------------------------------------
G4VHCIOentry* G4HCIOcatalog::GetEntry(const G4String& name)
{
  if(theCatalog.find(name) == theCatalog.end())
  {
    G4cout << "Hit Collection I/O manager entry \"" << name << "\" not found!"
           << std::endl;
    return nullptr;
  }
  else
  {
    G4VHCIOentry* ds = theCatalog[name];
    return ds;
  }
}

// --------------------------------------------------------------------
G4VPHitsCollectionIO* G4HCIOcatalog::GetHCIOmanager(const G4String& name)
{
  if(theStore.find(name) == theStore.cend())
  {
    G4cout << "Hit Collection I/O manager \"" << name << "\" not found!"
           << G4endl;
    return nullptr;
  }
  else
  {
    G4VPHitsCollectionIO* ds = theStore[name];
    return ds;
  }
}

// --------------------------------------------------------------------
void G4HCIOcatalog::PrintEntries()
{
  G4cout << "I/O manager entries: ";
  G4cout << theCatalog.size() << G4endl;
  for(auto it = theCatalog.cbegin(); it != theCatalog.cend(); ++it)
  {
    G4cout << "  --- " << (*it).first << G4endl;
  }
}

// --------------------------------------------------------------------
G4String G4HCIOcatalog::CurrentHCIOmanager()
{
  G4String list = "";
  for(auto it = theStore.cbegin(); it != theStore.cend(); ++it)
  {
    list += (*it).first + " ";
  }
  return list;
}

// --------------------------------------------------------------------
void G4HCIOcatalog::PrintHCIOmanager()
{
  G4cout << "I/O managers: ";
  G4cout << theStore.size() << G4endl;
  for(auto it = theStore.cbegin(); it != theStore.cend(); ++it)
  {
    G4cout << "  --- " << (*it).first << ", " << (*it).second << "." << G4endl;
  }
}

// --------------------------------------------------------------------
G4VPHitsCollectionIO* G4HCIOcatalog::GetHCIOmanager(G4int n)
{
  G4int i = 0;
  for(auto it = theStore.cbegin(); it != theStore.cend(); ++it)
  {
    if(i++ == n)
      return (*it).second;
  }
  return nullptr;
}
