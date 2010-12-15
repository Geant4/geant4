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
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL  and INFN Roma, gruppo collegato SanitÃ , Italy
// *Istituto Superiore di SanitÃ  and INFN Roma, gruppo collegato SanitÃ , Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2Main.hh"

#include "ML2PrimaryGenerationAction.hh"
#include "ML2WorldConstruction.hh"
#include "ML2PhysicsList.hh"
#include "ML2SteppingAction.hh"
#include "ML2EventAction.hh"
#include "ML2TrackingAction.hh"

#include "G4UImanager.hh"
#include "G4Timer.hh"
	

#ifdef G4VIS_USE
	#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
	#include "G4UIExecutive.hh"
#endif



#ifdef G4VIS_USE
void visio(int argc, char* argv[])
{
		G4VisManager *visManager=new G4VisExecutive;
		visManager->Initialize();
#ifdef G4UI_USE
		G4UIExecutive *ui = new G4UIExecutive(argc, argv);
                G4UImanager *UImanager = G4UImanager::GetUIpointer();
		UImanager->ApplyCommand("/control/execute vis.mac");     
		ui->SessionStart();
		delete ui;
#endif
		delete visManager;
}
#endif

int main(int argc, char* argv[])
{
	// instantiate the world class
	CML2WorldConstruction *myWorld=CML2WorldConstruction::GetInstance(); 

	G4RunManager *runManager=new G4RunManager();
	ML2PhysicsList *physics=new ML2PhysicsList();
	runManager->SetUserInitialization(physics);
	
	// build the primary generator
	CML2PrimaryGenerationAction *gun;
	gun = CML2PrimaryGenerationAction::GetInstance();

	// build the main messenger class for the input data
	CML2CInputData *myInputData;
	myInputData=new CML2CInputData();

	// initialize the primary generator variables
	gun->inizialize(&myInputData->inputData.primaryParticleData);

	// according to the number of the launching line
	if (argc==1)
	{
		myInputData->inputData.generalData.seed=1;
		myInputData->inputData.generalData.StartFileInputData="ml2.mac";
	}
	if (argc==2)
	{
		myInputData->inputData.generalData.seed=1;
		myInputData->inputData.generalData.StartFileInputData=(G4String)argv[1];
	}
	if (argc==3)
	{
		sscanf(argv[2],"%d", &myInputData->inputData.generalData.seed);
		myInputData->inputData.generalData.StartFileInputData=(G4String)argv[1];
	}

	// read the main mac file and execute the commands
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	G4String command = "/control/execute ";
	UImanager->ApplyCommand(command+myInputData->inputData.generalData.StartFileInputData); 


	// set and initialize the random generator
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	CLHEP::HepRandom :: setTheSeed(myInputData->inputData.generalData.seed);

	// create the world class
	if (!myWorld->create(&myInputData->inputData, myInputData->getbOnlyVisio()))
	{
		return 1; // if it fails to create the world
	}
	else
	{
		// initialize the primary generator according to the choosen particle source
		gun->design(myWorld->getCML2AcceleratorConstruction()->getAcceleratorIsoCentre());
	}

	// instantiate the convergence control class and assign it to myStepAction
	CML2Convergence *convergence=new CML2Convergence(myInputData->inputData.generalData.seed, myInputData->inputData.generalData.saving_in_Selected_Voxels_every_events, myInputData->inputData.generalData.fileExperimentalData, myInputData->inputData.generalData.fileExperimentalDataOut, myInputData->inputData.generalData.bCompareExp, myInputData->inputData.generalData.maxNumberOfEvents, gun->getNrecycling(), myInputData->inputData.generalData.nMaxLoop);


	// build the ML2RunAction to assign the single phantom name at each run 
	CML2RunAction *myRunAction=new CML2RunAction(convergence, myInputData->inputData.generalData.nBeam, myInputData->bOnlyVisio);
	
	CML2SteppingAction *myStepAction=new CML2SteppingAction(convergence);
	CML2EventAction *ML2EventAction = new CML2EventAction();

	runManager->SetUserInitialization(myWorld);
	runManager->SetUserAction(myRunAction);
	runManager->SetUserAction(gun);
	runManager->SetUserAction(myStepAction);
	runManager->SetUserAction(ML2EventAction);
	runManager->SetUserAction(new CML2TrackingAction);
	runManager->Initialize();

	// performances info 
	int nLoop=0;
	G4Timer MyFullTime;
	G4double loopElapsedTime;
	G4bool bStopRun=false;
	G4bool bNewGeometry=true;
	if (myInputData->bOnlyVisio)
	{
#ifdef G4VIS_USE
		// visualization
		myWorld->newGeometry();
		convergence->setNewGeometry();
		visio(argc, argv);
#endif
	}
	else
	{
		MyFullTime.Start();
		// compute 
		while (bNewGeometry)
		{
			bNewGeometry=myWorld->newGeometry();
			convergence->setNewGeometry();
			if (bNewGeometry)
			{
				if (CML2PhantomConstruction::GetInstance()->getPhantomName()!="Dicom1")
				{CML2WorldConstruction::GetInstance()->checkVolumeOverlap();}
				std::cout<<"################ START NEW GEOMETRY ########################"<<'\n';
				myRunAction->setActualLoop(nLoop);
				while (!bStopRun)
				{
					runManager->BeamOn(myInputData->inputData.generalData.nBeam);
					// check if the run has to be reapeted
					bStopRun=convergence->stopRun(); 
				}
				nLoop=0;
			std::cout<<"################ END NEW GEOMETRY ########################"<<'\n';
			}
			bStopRun=false;
		}
		MyFullTime.Stop();
		loopElapsedTime=MyFullTime.GetUserElapsed();
		std::cout << "loop elapsed time [s] : "<< loopElapsedTime << '\n';
		std::cout <<'\n';
	}
}

