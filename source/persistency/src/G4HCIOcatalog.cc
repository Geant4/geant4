// $Id: G4HCIOcatalog.cc,v 1.1 2002-11-24 13:45:24 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4HCIOcatalog.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4HCIOcatalog.hh"

// Addtional Include:
#include "G4VHCIOentry.hh"

G4HCIOcatalog* G4HCIOcatalog::f_thePointer = 0;

// Implementation of Constructor #1
G4HCIOcatalog::G4HCIOcatalog()
 : m_verbose(0)
{}

// Implementation of GetG4HCIOcatalog
G4HCIOcatalog* G4HCIOcatalog::GetG4HCIOcatalog()
{
  if ( f_thePointer == 0 ) f_thePointer = new G4HCIOcatalog;
  return f_thePointer;
}

// Implementation of RegisterEntry
void G4HCIOcatalog::RegisterEntry(G4VHCIOentry* d)
{
  if ( m_verbose > 0 ) {
    std::cout << "registering I/O manager entry \"" << d->GetName()
              << "\" " << d << "." << std::endl;
  }
  if ( theCatalog.find(d->GetName()) != theCatalog.end() ) {
    std::cout << "Redefining I/O Managers list " << d->GetName() << std::endl;
  } else {
    theCatalog[d->GetName()] = d;
  }
}

// Implementation of RegisterHCIOmanager
void G4HCIOcatalog::RegisterHCIOmanager(G4VPHitsCollectionIO* d)
{
  if ( m_verbose > 0 ) {
    std::cout << "registering I/O manager \"" << d->SDname()
              << "\" " << d << "." << std::endl;
  }
  if ( theStore.find(d->SDname()) != theStore.end() ) {
    std::cout << "Redefining I/O Manager " << d->SDname() << std::endl;
  } else {
    theStore[d->SDname()] = d;
  }
}

// Implementation of GetEntry
G4VHCIOentry* G4HCIOcatalog::GetEntry(std::string name)
{
  if ( theCatalog.find(name) == theCatalog.end() ) {
    std::cout << "Hit Collection I/O manager entry \"" << name
              << "\" not found!" << std::endl;
    return 0;
  } else {
    G4VHCIOentry* ds = theCatalog[name];
    return ds;
  }
}

// Implementation of GetHCIOmanager
G4VPHitsCollectionIO* G4HCIOcatalog::GetHCIOmanager(std::string name)
{
  if ( theStore.find(name) == theStore.end() ) {
    std::cout << "Hit Collection I/O manager \"" << name
              << "\" not found!" << std::endl;
    return 0;
  } else {
    G4VPHitsCollectionIO* ds = theStore[name];
    return ds;
  }
}

// Implementation of PrintEntries
void G4HCIOcatalog::PrintEntries()
{
  std::cout << "I/O manager entries: ";
  std::cout << theCatalog.size() << std::endl;
  HCIOmap::const_iterator it;
  for ( it=theCatalog.begin(); it != theCatalog.end(); it++ ) {
    std::cout << "  --- " << (*it).first << std::endl;
  }
}

// Implementation of CurrentHCIOmanager
std::string G4HCIOcatalog::CurrentHCIOmanager()
{
  std::string list = "";
  HCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    list += (*it).first + " ";
  }
  return list;
}

// Implementation of PrintHCIOmanager
void G4HCIOcatalog::PrintHCIOmanager()
{
  std::cout << "I/O managers: ";
  std::cout << theStore.size() << std::endl;
  HCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    std::cout << "  --- " << (*it).first
              << ", " << (*it).second << "." << std::endl;
  }
}

// Implementation of GetHCIOmanager
G4VPHitsCollectionIO* G4HCIOcatalog::GetHCIOmanager(int n)
{
  int i = 0;
  HCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    if (i++ == n) return (*it).second;
  }
  return 0;
}

// End of G4HCIOcatalog.cc

