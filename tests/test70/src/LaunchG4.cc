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
//      LaunchG4.cc
//
//      Copyright 2010 M. Karamitros <kara@cenbg.in2p3.fr>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include "LaunchG4.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "TrackingAction.hh"
#include "StackingAction.hh"
#include "ReactionAction.hh"
#include "G4AllITManager.hh"
#include "G4ITStepManager.hh"
#include "G4DNAChemistryManager.hh"

//#ifdef G4VIS_USE
//#include "G4VisExecutive.hh"
//#endif

#include "G4UImanager.hh"

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


using namespace std;

LaunchG4::LaunchG4()
{
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);

    fpRunManager = 0;
    fpPrimGenAct = 0;

#ifdef G4VIS_USE
//    fpVisManager = 0;
#endif

#ifdef G4UI_USE
    fpSession = 0;
#endif
}

LaunchG4::~LaunchG4()
{
#ifdef G4VIS_USE
//    if(fpVisManager)
//        delete fpVisManager;
#endif

#ifdef G4UI_USE
    if(fpSession)
        delete fpSession;
#endif

    if(fpRunManager)
        delete fpRunManager;
}

G4ParticleDefinition* LaunchG4::GetParticleDefinition()
{
    return fpPrimGenAct->GetParticleDefinition() ;
}

void LaunchG4::Initialize(G4double incidentEnergy, G4bool chemistryFlag)
{
    // Construct the default run manager
    fpRunManager = new G4RunManager;
    //------------------------------------------------------------------
    // Set mandatory user initialization classes
    fpRunManager->SetUserInitialization(new DetectorConstruction);

    // PhysicsList
    PhysicsList* thePhysicsList = new PhysicsList();
    thePhysicsList->SetChemistryFlag(chemistryFlag);
    fpRunManager->SetUserInitialization(thePhysicsList);

    //------------------------------------------------------------------
    // Set mandatory user action classes
    PrimaryGeneratorAction* primGenAction = new PrimaryGeneratorAction();

    if(incidentEnergy>0)
        primGenAction->SetIncidentEnergy(incidentEnergy);
    fpPrimGenAct = primGenAction ;

    fpRunManager->SetUserAction(primGenAction);

    //------------------------------------------------------------------
    // Set optional user action classes
    fpRunManager->SetUserAction(new RunAction());
    fpRunManager->SetUserAction(new TrackingAction());
    fpRunManager->SetUserAction(new SteppingAction(primGenAction));
    fpRunManager->SetUserAction(new StackingAction(chemistryFlag));
    if(chemistryFlag)
    {
        G4ITStepManager::Instance()->SetUserAction(new ReactionAction());
        G4ITStepManager::Instance()->SetVerbose(1);
    }

#ifdef G4VIS_USE
//    fpVisManager = new G4VisExecutive;
//    fpVisManager->Initialize();
#endif

    // Initialize G4 kernel
    fpRunManager->Initialize();

    G4DNAChemistryManager::Instance()->WriteInto("output.txt");
}

void LaunchG4::StartTimer()
{
//    if( clock_gettime( CLOCK_MONOTONIC, &fStart) == -1 )
//    {
//        perror( "clock gettime" );
//        exit(EXIT_FAILURE);
//    }
}

void LaunchG4::EndTimer()
{
//    if( clock_gettime( CLOCK_MONOTONIC, &fEnd) == -1 )
//    {
//        perror( "clock gettime" );
//        exit(EXIT_FAILURE);
//    }
}

void LaunchG4::RunSimu(G4bool drawing)
{
    // Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
    if(drawing)
    {
        G4cout << "Will start to process the script commands.mac" << G4endl;
        UI->ApplyCommand("/control/execute vis.mac");
    }
#endif

    G4cout << "Will start to process the script commands.mac" << G4endl;

    StartTimer();
    UI->ApplyCommand("/control/execute commands.mac");
    EndTimer();

    G4cout << "Has finished to process the script commands.mac" << G4endl;

#ifdef G4VIS_USE
//    if(drawing)
//        UI->ApplyCommand("/control/execute vis2.mac");
#endif


//    G4double executionTime = double(fEnd.tv_sec - fStart.tv_sec) + (fEnd.tv_nsec - fStart.tv_nsec)/1e9;
//    G4cout<<"Execution time =" << executionTime <<G4endl;
}

void LaunchG4::RunSimu(G4String macFile, G4bool /*drawing*/)
{
    // Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();

    G4String command = "/control/execute " ;
    command += macFile ;

    StartTimer();
    UI->ApplyCommand(command);
    EndTimer();

//    G4double executionTime = double(fEnd.tv_sec - fStart.tv_sec) + (fEnd.tv_nsec - fStart.tv_nsec)/1e9;
//    G4cout<<"Execution time =" << executionTime <<G4endl;
}

void LaunchG4::NewSession(int argc,char** argv)
{
#ifdef G4UI_USE
    // Define UI fpSession for interactive mode.
    // fpSession = new G4UIterminal(new G4UItcsh);
    if(fpSession == 0)
    {
#if defined(G4UI_USE_QT)
        G4cout << "You use Qt, please wait while processing,"
               << " if you chose the rich trajectory, this may take a long while, "
               << " but at the end, a window will appear ..." << G4endl;
#endif
        fpSession = new G4UIExecutive(argc,argv);
    }

    G4UImanager* UI = G4UImanager::GetUIpointer();
        if (fpSession->IsGUI())
          UI->ApplyCommand("/control/execute gui.mac");
#endif
}


void LaunchG4::StartSession()
{
#ifdef G4UI_USE
    // Define UI fpSession for interactive mode.
    if(fpSession)
    fpSession->SessionStart();
#endif
}
