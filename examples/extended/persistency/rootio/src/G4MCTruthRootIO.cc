// $Id: G4MCTruthRootIO.cc,v 1.2 2002-12-04 14:12:26 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4MCTruthRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4MCTruthRootIO.hh"

G4MCTruthRootIO* G4MCTruthRootIO::thePointer=G4MCTruthRootIO::GetMCTruthRootIO();

// Implementation of Store
bool G4MCTruthRootIO::Store(G4MCTEvent* mctevent)
{
  return true;
}

// Implementation of Retrieve
bool G4MCTruthRootIO::Retrieve(G4MCTEvent* & mctevent)
{
  return true;
}

// Implementation of GetMCTruthRootIO
G4MCTruthRootIO* G4MCTruthRootIO::GetMCTruthRootIO()
{
  return 0;
}

// End of G4MCTruthRootIO.cc

