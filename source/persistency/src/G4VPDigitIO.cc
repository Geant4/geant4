// $Id: G4VPDigitIO.cc,v 1.1 2002-11-24 13:45:25 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPDigitIO.cc
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#include "G4VPDigitIO.hh"

G4VPDigitIO* G4VPDigitIO::f_G4VPDigitIO = 0;

// Implementation of Constructor #1
G4VPDigitIO::G4VPDigitIO()
 : m_verbose(0)
{
  f_catalog = G4DCIOcatalog::GetG4DCIOcatalog();
}

// Implementation of SetVerboseLevel
void G4VPDigitIO::SetVerboseLevel(int v)
{
  m_verbose = v;

  // Loop through the registered Digit I/O managers
  for ( size_t i=0; i < f_catalog->NumberOfDCIOmanager(); i++ ) {
    G4VPDigitsCollectionIO* digitIOman = f_catalog->GetDCIOmanager(i);
    digitIOman->SetVerboseLevel(v);
  }
}

// End of G4VPDigitIO.cc

