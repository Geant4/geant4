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
#include "ML2Main.hh"

#include "ML2PrimaryGenerationAction.hh"
#include "ML2WorldConstruction.hh"
#include "ML2PhysicsList.hh"
#include "ML2SteppingAction.hh"
#include "ML2EventAction.hh"
#include "ML2TrackingAction.hh"

#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "G4Timer.hh"
	
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


void visio(int argc, char* argv[])
{
	G4VisManager *visManager=new G4VisExecutive;
	visManager -> Initialize();
	G4UIExecutive *ui = new G4UIExecutive(argc, argv);
			G4UImanager *UImanager = G4UImanager::GetUIpointer();
	UImanager -> ApplyCommand("/control/execute macroAndData/macro_files/main/vis.mac");
	ui -> SessionStart();
	delete ui;
	delete visManager;
}

int main(int argc, char* argv[])
{
	G4RunManager *runManager = new G4RunManager();
	// instantiate the world class
	CML2WorldConstruction *myWorld=CML2WorldConstruction::GetInstance(); 

	// read the main mac file and execute the commands
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
        UImanager->ApplyCommand("/control/verbose 2");
        UImanager->ApplyCommand("/run/verbose 2");

	ML2PhysicsList *physics = new ML2PhysicsList();
	runManager -> SetUserInitialization(physics);

	// instantiate the scoring manager
	G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
	scoringManager->SetVerboseLevel(1);
	
	// build the primary generator
	CML2PrimaryGenerationAction *gun;
	gun = CML2PrimaryGenerationAction::GetInstance();

	// build the main messenger class for the input data
	CML2CInputData *myInputData = new CML2CInputData();

	// initialize the primary generator variables
	gun->inizialize(&myInputData->inputData.primaryParticleData);

	// according to the number of the launching line
	if (argc==1)
	{
	  myInputData->inputData.generalData.seed = 1;
	  myInputData->inputData.generalData.StartFileInputData = "macroAndData/macro_files/main/ml2.mac";
	}
	if (argc==2)
	{
	  myInputData->inputData.generalData.seed = 1;
	  myInputData->inputData.generalData.StartFileInputData = (G4String)argv[1];
	}
	if (argc==3)
	{
	  sscanf(argv[2],"%d", &myInputData->inputData.generalData.seed);
	  myInputData->inputData.generalData.StartFileInputData = (G4String)argv[1];
	}

	G4String command = "/control/execute ";
	UImanager->ApplyCommand(command+myInputData->inputData.generalData.StartFileInputData);


	// set and initialize the random generator
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	CLHEP::HepRandom::setTheSeed(myInputData->inputData.generalData.seed);
	G4cout << "Using seed " << CLHEP::HepRandom::getTheSeed() << G4endl;

	// create the world class
	if (!myWorld->create(&myInputData->inputData, myInputData->getbOnlyVisio()))
	{
	  return 1; // if it fails to create the world
	}
	else
	{
	  // initialize the primary generator according to the chosen particle source
	  gun->design(myWorld->getCML2AcceleratorConstruction()->getAcceleratorIsoCentre());
        }

	// instantiate the convergence control class and assign it to myStepAction
	CML2Convergence *convergence = new CML2Convergence(myInputData->inputData.generalData.seed,
			myInputData->inputData.generalData.saving_in_Selected_Voxels_every_events,
			myInputData->inputData.generalData.fileExperimentalData,
			myInputData->inputData.generalData.fileExperimentalDataOut,
			myInputData->inputData.generalData.bCompareExp,
			myInputData->inputData.generalData.maxNumberOfEvents,
			gun->getNrecycling(),
			myInputData->inputData.generalData.nMaxLoop);


	// build the ML2RunAction to assign the single phantom name at each run 
	CML2RunAction *myRunAction = new CML2RunAction(convergence, myInputData->inputData.generalData.nBeam,
			myInputData->bOnlyVisio,
                        myInputData->inputData.voxelSegmentation.nX,
                        myInputData->inputData.voxelSegmentation.nY,
                        myInputData->inputData.voxelSegmentation.nZ );
	
	CML2SteppingAction *myStepAction = new CML2SteppingAction(convergence);
	CML2EventAction *ML2EventAction  = new CML2EventAction();

	runManager -> SetUserInitialization(myWorld);
	runManager -> SetUserAction(myRunAction);
	runManager -> SetUserAction(gun);
	runManager -> SetUserAction(myStepAction);
	runManager -> SetUserAction(ML2EventAction);
	runManager -> SetUserAction(new CML2TrackingAction);
	runManager -> Initialize();

	// performances info 
	int nLoop = 0;
	G4Timer MyFullTime;
	G4double loopElapsedTime;
	G4bool bStopRun = false;
	G4bool bNewGeometry = true;

	if (myInputData->bOnlyVisio)
	{
	  // visualization
	  myWorld -> newGeometry();
	  convergence -> setNewGeometry();
	  visio(argc, argv);
	}
	else
	{
	  MyFullTime.Start();
	  // compute 
	  while (bNewGeometry)
	    {
	      bNewGeometry = myWorld -> newGeometry();
	      convergence -> setNewGeometry();

	      if (bNewGeometry)
		{
		  runManager -> Initialize();
		  CML2WorldConstruction::GetInstance()->checkVolumeOverlap();

		  G4cout << "################ START NEW GEOMETRY ########################" << G4endl;
		  myRunAction->setActualLoop(nLoop);

		  while (!bStopRun)
		    {
		      runManager->BeamOn(myInputData->inputData.generalData.nBeam);
		      // check if the run has to be repeated
		      bStopRun = convergence->stopRun();
		    }
		  nLoop = 0;
		  G4cout << "################ END NEW GEOMETRY ########################" << G4endl;
		}
	      bStopRun = false;
	    }
	  MyFullTime.Stop();
	  loopElapsedTime = MyFullTime.GetUserElapsed();
	  G4cout << "loop elapsed time [s] : "<< loopElapsedTime << '\n' << G4endl;
	}
        delete runManager;
        return 0;
}

