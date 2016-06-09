//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RunAction.cc,v 1.2 2007-12-10 16:28:54 gunter Exp $
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

