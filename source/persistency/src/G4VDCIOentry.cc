// $Id: G4VDCIOentry.cc,v 1.2 2002-12-04 10:25:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VDCIOentry.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4VDCIOentry.hh"

// Addtional Include:
#include "G4DCIOcatalog.hh"

// Implementation of Constructor #1
G4VDCIOentry::G4VDCIOentry(G4std::string n)
 : m_name(n)
{
  G4DCIOcatalog* c = G4DCIOcatalog::GetG4DCIOcatalog();
  c->RegisterEntry(this);

  m_verbose = G4PersistencyCenter::GetG4PersistencyCenter()->VerboseLevel();
}

// End of G4VDCIOentry.cc

