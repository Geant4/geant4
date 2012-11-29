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
//    *******************************
//    *                             *
//    *    RemSimRunAction.cc       *
//    *                             *
//    *******************************
//
// Code developed by: S.Guatelli, susanna@uow.edu.au
// $Id$
//

#include "RemSimRunAction.hh"
#include "RemSimDetectorConstruction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "RemSimRunAction.hh"
#include "RemSimAnalysis.hh"

RemSimRunAction::RemSimRunAction()
{
 // Create analysis manager
 G4AnalysisManager::Instance();
}

RemSimRunAction::~RemSimRunAction()
{
 delete G4AnalysisManager::Instance();
 }

void RemSimRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun -> GetRunID() << " start." << G4endl;
 
 // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Open a ROOT output file
  //
  G4String fileName = "remsim.root";
  analysisManager -> OpenFile(fileName);
  analysisManager -> SetFirstHistoId(1);
  
  // Create the histograms
  analysisManager -> CreateH1("1","Energy of secondary p reaching the phantom",1000, 0., 10000.); 
  analysisManager -> CreateH1("2","Energy of secondary n reaching the phantom", 1000, 0., 10000.);  
  analysisManager -> CreateH1("3","Energy of secondary pions reaching the phantom",1000, 0., 10000.); 
  analysisManager -> CreateH1("4","Energy of secondary alpha reaching the phantom",100, 0., 100.); 
  analysisManager -> CreateH1("5","Energy of secondary p produced in the phantom",100, 0., 1000.); 
  analysisManager -> CreateH1("6","Energy of secondary n produced in the phantom", 100, 0., 1000.); 
  analysisManager -> CreateH1("7","Energy of secondary pions produced in the phantom",200, 0., 2000.);
  analysisManager -> CreateH1("8","Energy of secondary alpha produced in the phantom", 100, 0.,100.);  
}

void RemSimRunAction::EndOfRunAction(const G4Run* aRun)
{  
 G4double numberEvents = aRun -> GetNumberOfEvent();
 G4cout<< "Number of events:" << numberEvents << G4endl;

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Save the histograms and close the ROOT file
  analysisManager -> Write();
  analysisManager -> CloseFile();
}

