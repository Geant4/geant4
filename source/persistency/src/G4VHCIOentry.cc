// $Id: G4VHCIOentry.cc,v 1.3 2002-12-04 13:57:30 morita Exp $
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
G4VHCIOentry::G4VHCIOentry(G4std::string n)
 : m_name(n)
{
  G4HCIOcatalog* c = G4HCIOcatalog::GetHCIOcatalog();
  c->RegisterEntry(this);

  m_verbose = G4PersistencyCenter::GetPersistencyCenter()->VerboseLevel();
}

// End of G4VHCIOentry.cc

