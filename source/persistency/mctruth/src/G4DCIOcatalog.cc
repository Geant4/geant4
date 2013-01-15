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
// File: G4DCIOcatalog.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4DCIOcatalog.hh"

// Addtional Include:
#include "G4VDCIOentry.hh"

G4ThreadLocal G4DCIOcatalog* G4DCIOcatalog::f_thePointer = 0;

// Implementation of Constructor #1
G4DCIOcatalog::G4DCIOcatalog()
 : m_verbose(0)
{}

// Implementation of GetDCIOcatalog
G4DCIOcatalog* G4DCIOcatalog::GetDCIOcatalog()
{
  if ( f_thePointer == 0 ) f_thePointer = new G4DCIOcatalog;
  return f_thePointer;
}

// Implementation of RegisterEntry
void G4DCIOcatalog::RegisterEntry(G4VDCIOentry* d)
{
  if ( m_verbose > 0 ) {
    G4cout << "registering I/O manager entry \"" << d->GetName()
              << "\" " << d << "." << G4endl;
  }
  if ( theCatalog.find(d->GetName()) != theCatalog.end() ) {
    G4cout << "Redefining I/O Managers list " << d->GetName() << G4endl;
  } else {
    theCatalog[d->GetName()] = d;
  }
}

// Implementation of RegisterDCIOmanager
void G4DCIOcatalog::RegisterDCIOmanager(G4VPDigitsCollectionIO* d)
{
  if ( m_verbose > 0 ) {
    G4cout << "registering I/O manager \"" << d->DMname()
              << "\" " << d << "." << G4endl;
  }
  if ( theStore.find(d->DMname()) != theStore.end() ) {
    G4cout << "Redefining I/O Manager " << d->DMname() << G4endl;
  } else {
    theStore[d->DMname()] = d;
  }
}

// Implementation of GetEntry
G4VDCIOentry* G4DCIOcatalog::GetEntry(std::string name)
{
  if ( theCatalog.find(name) == theCatalog.end() ) {
    G4cout << "Digit Collection I/O manager entry \"" << name
              << "\" not found!" << G4endl;
    return 0;
  } else {
    G4VDCIOentry* ds = theCatalog[name];
    return ds;
  }
}

// Implementation of GetDCIOmanager
G4VPDigitsCollectionIO* G4DCIOcatalog::GetDCIOmanager(std::string name)
{
  if ( theStore.find(name) == theStore.end() ) {
    G4cout << "Digit Collection I/O manager \"" << name
              << "\" not found!" << G4endl;
    return 0;
  } else {
    G4VPDigitsCollectionIO* ds = theStore[name];
    return ds;
  }
}

// Implementation of PrintEntries
void G4DCIOcatalog::PrintEntries()
{
  G4cout << "I/O manager entries: ";
  G4cout << theCatalog.size() << G4endl;
  DCIOmap::const_iterator it;
  for ( it=theCatalog.begin(); it != theCatalog.end(); it++ ) {
    G4cout << "  --- " << (*it).first << G4endl;
  }
}

// Implementation of CurrentDCIOmanager
std::string G4DCIOcatalog::CurrentDCIOmanager()
{
  std::string list = "";
  DCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    list += (*it).first + " ";
  }
  return list;
}

// Implementation of PrintDCIOmanager
void G4DCIOcatalog::PrintDCIOmanager()
{
  G4cout << "I/O managers: ";
  G4cout << theStore.size() << G4endl;
  DCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    G4cout << "  --- " << (*it).first
           << ", " << (*it).second << "." << G4endl;
  }
}

// Implementation of GetDCIOmanager
G4VPDigitsCollectionIO* G4DCIOcatalog::GetDCIOmanager(int n)
{
  int i = 0;
  DCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    if (i++ == n) return (*it).second;
  }
  return 0;
}

// End of G4DCIOcatalog.cc

