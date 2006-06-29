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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************//

#include "globals.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "DicomPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "DicomGeometry.hh"
#include "DicomPrimaryGeneratorAction.hh"
#include "DicomEventAction.hh"
#include "DicomHandler.hh"

int main(int argc,char** argv)
{
				
  // Treatment of DICOM images before creating the G4runManager
  DicomHandler* dcmHandler = new DicomHandler;
  dcmHandler->checkFileFormat();

  // Initialisation of physics, geometry, primary particles ... 
  G4RunManager* runManager = new G4RunManager;
  DicomGeometry* theGeometry = new DicomGeometry();
  runManager->SetUserInitialization(new DicomPhysicsList);
  runManager->SetUserInitialization(theGeometry);
  runManager->SetUserAction(new DicomPrimaryGeneratorAction());
  runManager->SetUserAction(new DicomEventAction);

  runManager->Initialize();

#ifdef G4VIS_USE
  // visualisation manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  G4UImanager* UI = G4UImanager::GetUIpointer();
 
  if (argc==1)
    {
      G4UIsession* session = new G4UIterminal(new G4UItcsh);
      UI->ApplyCommand("/control/execute default.mac");
      session->SessionStart();
      delete session;
    }
  else
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  delete runManager;

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete dcmHandler;

  return 0;
}




