// $Id: G4DigitRootIO.cc,v 1.2 2002-12-04 14:12:26 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4DigitRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4DigitRootIO.hh"

// Implementation of Constructor #1
G4DigitRootIO::G4DigitRootIO()
{
  f_G4VPDigitIO = (G4VPDigitIO*) this;
}

// Implementation of Destructor #1
G4DigitRootIO::~G4DigitRootIO()
{}

// Implementation of GetDigitRootIO
G4DigitRootIO* G4DigitRootIO::GetDigitRootIO()
{
  G4DigitRootIO* dio=0;
  if ( f_G4VPDigitIO == 0 ) {
    dio = new G4DigitRootIO;
  }
  return dio;
}

// Implementation of Store
bool G4DigitRootIO::Store(const G4DCofThisEvent* dcevt)
{
  return true;
}

// Implementation of Retrieve
bool G4DigitRootIO::Retrieve(G4DCofThisEvent*& dcevt)
{
  return true;
}

// End of G4DigitRootIO.cc

