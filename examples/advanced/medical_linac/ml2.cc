
#include "ML2Main.h"

#include "ML2PrimaryGenerationAction.h"
#include "ML2WorldConstruction.h"
#include "ML2PhysicsList.h"
#include "ML2SteppingAction.h"
#include "ML2EventAction.h"
#include "ML2TrackingAction.h"

#include "G4UImanager.hh"
#include "G4Timer.hh"
	

#include "G4UIterminal.hh"
#ifdef G4UI_USE_XM
	#include "G4UIXm.hh"
#endif
#ifdef G4VIS_USE
	#include "G4VisExecutive.hh"
#endif



void visio()
{
	G4UIsession *session=0;
	session=new G4UIterminal();
	#ifdef G4VIS_USE
		G4VisManager * visManager=new G4VisExecutive;
		visManager->Initialize();
	#endif
		G4UImanager *UI=G4UImanager::GetUIpointer();
		if (session)
		{
			UI->ApplyCommand("/control/execute vis.mac");
			session->SessionStart();
			delete session;
		}

#ifdef G4VIS_USE
 delete visManager;
#endif
}


int main(int argc, char* argv[])
{
	if (argc==3)
	{
// record launch data
		CML2CInputData *myInputData;
		myInputData=new CML2CInputData();
		sscanf(argv[2],"%d", &myInputData->inputData.generalData.seed);
		myInputData->inputData.generalData.StartFileInputData=(G4String)argv[1];

		G4RunManager *runManager=new G4RunManager();
		CML2PhysicsList *physics=new CML2PhysicsList();
		runManager->SetUserInitialization(physics);
		CML2MainMessenger *ML2MainMessenger;
		ML2MainMessenger=new CML2MainMessenger(myInputData);

		CML2WorldConstruction *myWorld=CML2WorldConstruction::GetInstance(); 
	
		CML2PrimaryGenerationAction *gun;
		gun = new CML2PrimaryGenerationAction(&myInputData->inputData.primaryParticleData);

		G4UImanager* UI = G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		UI->ApplyCommand(command+myInputData->inputData.generalData.StartFileInputData); 

// initialize the primary generator according to the choosen particle source
		gun->design();

		CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
		CLHEP::HepRandom :: setTheSeed(myInputData->inputData.generalData.seed);

		myWorld->create(&myInputData->inputData);
		
		runManager->SetUserInitialization(myWorld);


		runManager->SetUserAction(gun);

// instantiate the convergence control class and assign it to myStepAction
		CML2Convergence *convergence=new CML2Convergence(myInputData->inputData.generalData.seed, myInputData->inputData.generalData.saving_in_Selected_Voxels_every_events, myInputData->inputData.generalData.fileExperimentalData, myInputData->inputData.generalData.bCompareExp, myInputData->inputData.generalData.minNumberOfEvents);
		CML2SteppingAction *myStepAction=new CML2SteppingAction(convergence);
		runManager->SetUserAction(myStepAction);

		CML2EventAction *ML2EventAction = new CML2EventAction();
		runManager->SetUserAction(ML2EventAction);

		runManager->SetUserAction(new CML2TrackingAction);

		runManager->Initialize();

		int nLoop=0;
		G4Timer MyTime, MyTimeStop;
		G4double loopElapsedTime;

		G4bool bAgain=true;
		if (myInputData->bOnlyVisio)
		{
			visio();
		}
		else
		{
			MyTimeStop.Start();
			while (bAgain)
			{
				std::cout<<"*********************************************"<<'\n';
				std::cout << "loop n. "<<++nLoop<<'\n';
					MyTime.Start();
					std::cout << "Launched "<< myInputData->inputData.generalData.nBeam <<" random primary particles" << '\n';

					runManager->BeamOn(myInputData->inputData.generalData.nBeam);

					MyTime.Stop();
					loopElapsedTime=MyTime.GetUserElapsed();
					std::cout << "loop elapsed time [s] : "<< loopElapsedTime << '\n';
					std::cout <<'\n';
// check if the run have to be reapeted
				bAgain=convergence->runAgain();
			}
		}
		//job termination
		delete myInputData;
//		delete G4MainMessenger;
		delete myWorld;
//		delete runManager;
	}
	return 0;
}

