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
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
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
	CML2CInputData *myInputData;
	myInputData=new CML2CInputData();
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

		G4RunManager *runManager=new G4RunManager();
		ML2PhysicsList *physics=new ML2PhysicsList();
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
	return 0;
}

