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
//
// $Id: Test2Run.cc,v 1.5 2010-11-03 08:48:57 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "Test2Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

Test2Run::Test2Run() {

  G4String detMassSD = "MassWorldSD";
  G4String detParaSD = "ParallelWorldSD";
  G4String hcname    = "PhantomCollection";
  //
  G4String detMassPS = "MassWorldPS";
  G4String detParaPS = "ParallelWorldPS";
  G4String primNameSum[16] = {"eDep",
			     "trackLengthGamma",
			     "trackLengthElec",
			     "trackLengthPosi",
			     "nStepGamma",
			     "nStepElec",
			     "nStepPosi",
			     "PassageTrackLength",
			     "DoseDeposit",
			     "FlatSurfaceCurrent",
			     "PassageCellCurrent",
			     "FlatSurfaceFlux",
			     "CellFlux",
			     "PassageCellFlux",
			     "NofSecondary",
			     "CellCharge" };

  std::vector<G4String> prim;
  for ( G4int i = 0; i < 16; i++){
    prim.push_back(primNameSum[i]);
  }

  MassRunSD = new Test2RunSD(detMassSD,hcname,prim);
  ParaRunSD = new Test2RunSD(detParaSD,hcname,prim);

  MassRunPS = new Test2RunPS(detMassPS,prim);
  ParaRunPS = new Test2RunPS(detParaPS,prim);

}

Test2Run::~Test2Run() {
  if ( MassRunSD ) delete MassRunSD;
  if ( ParaRunSD ) delete ParaRunSD;
  if ( MassRunPS ) delete MassRunPS;
  if ( ParaRunPS ) delete ParaRunPS;
}

#include "Test2PhantomHit.hh"

void Test2Run::RecordEvent(const G4Event* evt) {

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;

  MassRunSD->RecordEvent(evt);
  ParaRunSD->RecordEvent(evt);
  MassRunPS->RecordEvent(evt);
  ParaRunPS->RecordEvent(evt);
}

void Test2Run::DumpQuantitiesToFile() {
  MassRunSD->DumpQuantitiesToFile();
  ParaRunSD->DumpQuantitiesToFile();
  MassRunPS->DumpQuantitiesToFile();
  ParaRunPS->DumpQuantitiesToFile();
}
