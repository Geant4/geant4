// $Id: G4VPHitIO.cc,v 1.1 2002-11-24 13:45:25 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPHitIO.cc
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#include "G4VPHitIO.hh"

G4VPHitIO* G4VPHitIO::f_G4VPHitIO = 0;

// Implementation of Constructor #1
G4VPHitIO::G4VPHitIO()
 : m_verbose(0)
{
  f_catalog = G4HCIOcatalog::GetG4HCIOcatalog();
}

// Implementation of SetVerboseLevel
void G4VPHitIO::SetVerboseLevel(int v)
{
  m_verbose = v;

  // Loop through the registered Hit I/O managers
  for ( size_t i=0; i < f_catalog->NumberOfHCIOmanager(); i++ ) {
    G4VPHitsCollectionIO* hitIOman = f_catalog->GetHCIOmanager(i);
    hitIOman->SetVerboseLevel(v);
  }
}

// End of G4VPHitIO.cc

