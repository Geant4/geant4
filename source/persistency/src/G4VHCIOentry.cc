// $Id: G4VHCIOentry.cc,v 1.1 2002-11-24 13:45:24 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VHCIOentry.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4VHCIOentry.hh"

// Addtional Include:
#include "G4HCIOcatalog.hh"

// Implementation of Constructor #1
G4VHCIOentry::G4VHCIOentry(std::string n)
 : m_name(n)
{
  G4HCIOcatalog* c = G4HCIOcatalog::GetG4HCIOcatalog();
  c->RegisterEntry(this);

  m_verbose = G4PersistencyCenter::GetG4PersistencyCenter()->VerboseLevel();
}

// End of G4VHCIOentry.cc

