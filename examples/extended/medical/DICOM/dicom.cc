//*******************************************************
// G4 routine to build a geometry based on DICOM
// images for medical purpose.  Material are extracted
// from CT (Hounsfield) numbers and built using ICRU
// report 46.
//
// Note:-DicomGeometry (DetectorConstruction), DicomHandler,
//   	 DicomPatientParametrisation .cc and .hh files are
//	 the ones necessary for the correct construction of
//	 a geometry based on .dcm files
//	-The PrimaryGeneratorAction is designed to send
//	 irradiation in a way similar to medical irradiation
//	 but it can be replaced without any problems
//	-The PhysicsList is designed for the low energy
//
// WARNING : this version of the code is slow :
//	     1 minute for 100 initial photons, this will
//	     be corrected in the future
//
// The code was written by :
//	Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
// For more information contact :
//	Louis Archambault louis.archambault@phy.ulaval.ca
// or
//	Luc Beaulieu beaulieu@phy.ulaval.ca
//
// Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//*******************************************************

#include "globals.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "DicomPhysicsList.hh"
#include "DicomVisManager.hh"
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

  // visualisation manager
  G4VisManager* visManager = new DicomVisManager;
  visManager->Initialize();

  HepRandom::setTheEngine(new RanecuEngine);

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
  delete visManager;
  delete dcmHandler;

  return 0;
}




