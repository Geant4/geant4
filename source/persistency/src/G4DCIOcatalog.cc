// $Id: G4DCIOcatalog.cc,v 1.1 2002-11-24 13:45:24 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4DCIOcatalog.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4DCIOcatalog.hh"

// Addtional Include:
#include "G4VDCIOentry.hh"

G4DCIOcatalog* G4DCIOcatalog::f_thePointer = 0;

// Implementation of Constructor #1
G4DCIOcatalog::G4DCIOcatalog()
 : m_verbose(0)
{}

// Implementation of GetG4DCIOcatalog
G4DCIOcatalog* G4DCIOcatalog::GetG4DCIOcatalog()
{
  if ( f_thePointer == 0 ) f_thePointer = new G4DCIOcatalog;
  return f_thePointer;
}

// Implementation of RegisterEntry
void G4DCIOcatalog::RegisterEntry(G4VDCIOentry* d)
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

// Implementation of RegisterDCIOmanager
void G4DCIOcatalog::RegisterDCIOmanager(G4VPDigitsCollectionIO* d)
{
  if ( m_verbose > 0 ) {
    std::cout << "registering I/O manager \"" << d->DMname()
              << "\" " << d << "." << std::endl;
  }
  if ( theStore.find(d->DMname()) != theStore.end() ) {
    std::cout << "Redefining I/O Manager " << d->DMname() << std::endl;
  } else {
    theStore[d->DMname()] = d;
  }
}

// Implementation of GetEntry
G4VDCIOentry* G4DCIOcatalog::GetEntry(std::string name)
{
  if ( theCatalog.find(name) == theCatalog.end() ) {
    std::cout << "Digit Collection I/O manager entry \"" << name
              << "\" not found!" << std::endl;
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
    std::cout << "Digit Collection I/O manager \"" << name
              << "\" not found!" << std::endl;
    return 0;
  } else {
    G4VPDigitsCollectionIO* ds = theStore[name];
    return ds;
  }
}

// Implementation of PrintEntries
void G4DCIOcatalog::PrintEntries()
{
  std::cout << "I/O manager entries: ";
  std::cout << theCatalog.size() << std::endl;
  DCIOmap::const_iterator it;
  for ( it=theCatalog.begin(); it != theCatalog.end(); it++ ) {
    std::cout << "  --- " << (*it).first << std::endl;
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
  std::cout << "I/O managers: ";
  std::cout << theStore.size() << std::endl;
  DCIOstore::const_iterator it;
  for ( it=theStore.begin(); it != theStore.end(); it++ ) {
    std::cout << "  --- " << (*it).first
              << ", " << (*it).second << "." << std::endl;
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

