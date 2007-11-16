// $Id: RunAction.cc,v 1.1 2007-11-16 14:29:33 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   RunAction.cc
//
//                                         2007 Q
// ====================================================================
#include "RunAction.hh"
#include "Analysis.hh"
#include "G4MPImanager.hh"
#include <stdio.h>

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////
RunAction::RunAction()
//////////////////////////
{
}


///////////////////////////
RunAction::~RunAction()
///////////////////////////
{
}


//////////////////////////////////////////////
void RunAction::BeginOfRunAction(const G4Run*)
//////////////////////////////////////////////
{
  Analysis* myana= Analysis::GetAnalysis();
  myana-> Clear();
}


////////////////////////////////////////////
void RunAction::EndOfRunAction(const G4Run*)
////////////////////////////////////////////
{
  G4int rank= G4MPImanager::GetManager()-> GetRank();
  
  char str[64];
  sprintf(str, "dose-%03d.root", rank);
  G4String fname(str);

  Analysis* myana= Analysis::GetAnalysis();
  myana-> Save(fname);
}

