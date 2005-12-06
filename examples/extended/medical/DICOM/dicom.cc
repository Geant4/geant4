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




