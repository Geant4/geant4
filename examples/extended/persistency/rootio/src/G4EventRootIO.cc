// $Id: G4EventRootIO.cc,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4EventRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4EventRootIO.hh"

// Implementation of Constructor #1
G4EventRootIO::G4EventRootIO( G4PersistencyManager* pc )
 : f_pc(pc), f_hepmcio(0), f_mctruthio(0), f_PHCMan(0), f_PDCMan(0)
{
  f_hepmcio   = (G4HepMCRootIO*)   f_pc->HepMCIO();
  f_mctruthio = (G4MCTruthRootIO*) f_pc->MCTruthIO();
  f_PHCMan    = (G4HitRootIO*)     f_pc->HitIO();
  f_PDCMan    = (G4DigitRootIO*)   f_pc->DigitIO();
}

// Implementation of Store
// G4bool G4EventRootIO::Store( HepMC::GenEvent* hepevt, G4MCTEvent* mctevt, const G4Event* g4evt)
G4bool G4EventRootIO::Store( HepMC::GenEvent* hepevt, const G4Event* g4evt)
{
  // if (hepevt == 0) return false; // no event generated...

  // G4RootEvent* evt = new G4RootEvent();

  return true;
}

// Implementation of Store
G4bool G4EventRootIO::Store( const G4Event* anEvent )
{
  return true;
}

// Implementation of Retrieve
G4bool G4EventRootIO::Retrieve( G4Pevent*& anEvent )
{
  return true;
}

// Implementation of Retrieve
G4bool G4EventRootIO::Retrieve( G4Event*& anEvent )
{
  return true;
}

// End of G4EventRootIO.cc

