//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StandardPrintableScorer.cc,v 1.1 2002-07-29 16:03:15 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4StandardPrintableScorer.cc
//
// ----------------------------------------------------------------------
#include "G4StandardPrintableScorer.hh"

#include "G4StandardScorer.hh"
#include "G4StandardScoreTable.hh"

G4StandardPrintableScorer::G4StandardPrintableScorer():
  fStandardScorer(new G4StandardScorer)
{}
G4StandardPrintableScorer::~G4StandardPrintableScorer(){
  delete fStandardScorer;
}

G4VPScorer *G4StandardPrintableScorer::GetPointerToScorer(){
  return fStandardScorer;
}

void G4StandardPrintableScorer::
Score(const G4Step &step, const G4PStep &pstep){
  fStandardScorer->Score(step, pstep);
}

void G4StandardPrintableScorer::
PrintTable(G4std::ostream *out,
	   const G4VIStore *istore){
  G4StandardScoreTable sc_table(istore);
  sc_table.Print(fStandardScorer->GetMapPtkStandardCellScorer(), out);
}
