// $Id: G4MCTruthRootIO.cc,v 1.1 2002-12-04 02:44:29 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4MCTruthRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4MCTruthRootIO.hh"

G4MCTruthRootIO* G4MCTruthRootIO::thePointer=G4MCTruthRootIO::GetG4MCTruthRootIO();

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

// Implementation of GetG4MCTruthRootIO
G4MCTruthRootIO* G4MCTruthRootIO::GetG4MCTruthRootIO()
{
  return 0;
}

// End of G4MCTruthRootIO.cc

