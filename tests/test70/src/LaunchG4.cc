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
#include "SteppingVerbose.hh"
#include "ActionInitialization.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4DNAChemistryManager.hh"
#include <iomanip>
#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include "Parser.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

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
  fpVisManager = 0;
#endif

#ifdef G4UI_USE
  fpSession = 0;
#endif
}

LaunchG4::~LaunchG4()
{
#ifdef G4VIS_USE
  if(fpVisManager)
  delete fpVisManager;
#endif

#ifdef G4UI_USE
  if(fpSession)
  delete fpSession;
#endif

  if (fpRunManager) delete fpRunManager;
}

void LaunchG4::Initialize(G4bool chemistryFlag)
{

#ifdef G4VIS_USE
  fpVisManager = new G4VisExecutive;
  fpVisManager->Initialize();
#endif

  Command* commandLine = 0;
  if ((commandLine = CommandLineParser::GetParser()->GetCommandIfActive("-mt")))
  {

#ifdef G4MULTITHREADED

    fpRunManager= new G4MTRunManager;

    if(commandLine->fOption.empty())
    {
      ((G4MTRunManager*)fpRunManager)->SetNumberOfThreads(1);
    }
    else
    {
      int nThreads = G4UIcommand::ConvertToInt(commandLine->fOption);
      ((G4MTRunManager*)fpRunManager)->SetNumberOfThreads(nThreads);
    }
#else
    fpRunManager = new G4RunManager();

#endif

  }
  else
  {
    fpRunManager = new G4RunManager;
  }

  /**
   * Tells to the chemistry manager whether the chemistry
   * needs to be activated.
   * WARNING : if you don't use the chemistry do not activate it
   * otherwise it might generate memory leaks with tracks created but
   * not destroyed.
   */
  G4DNAChemistryManager::Instance()->SetChemistryActivation(chemistryFlag);

  //------------------------------------------------------------------
  // Set mandatory user initialization classes
  fpRunManager->SetUserInitialization(new DetectorConstruction);
  fpRunManager->SetUserInitialization(new PhysicsList());
  fpRunManager->SetUserInitialization(new ActionInitialization());

  G4DNAChemistryManager::Instance()->InitializeMaster();

  // Initialize G4 kernel
  fpRunManager->Initialize();
}

void LaunchG4::RunSimu(G4String macFile)
{
  // Get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  G4String command = "/control/execute ";
  command += macFile;

  G4cout << setfill('-') << setw(7) << " ";
  G4cout << "LaunchG4::RunSimu - Start script : " << macFile;
  G4cout << setfill('-') << setw(7) << " ";
  G4cout << G4endl;

  G4Timer timer;
  timer.Start();
  UI->ApplyCommand(command);
  timer.Stop();

  G4cout << setfill('-') << setw(7) << " ";
  G4cout << "LaunchG4::RunSimu - End script : " << macFile;
  G4cout << setfill('-') << setw(7) << " ";
  G4cout << G4endl;

  G4cout << setfill('-') << setw(7) << " ";
  G4cout << "User Elapsed Time (s) =" << timer.GetUserElapsed() / s << " ";
  G4cout << setfill('-') << setw(7) << " ";
  G4cout << G4endl;
}

void LaunchG4::NewSession(int argc, char** argv, const G4String& sessionType)
{
#ifdef G4UI_USE
  // Define UI fpSession for interactive mode.
  if(fpSession == 0)
  {
    fpSession = new G4UIExecutive(argc,argv,sessionType);
    BuildGUIFrame();
  }
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

void LaunchG4::BuildGUIFrame()
{
#ifdef G4UI_USE
  if(fpSession && fpSession -> IsGUI())
  {
    G4UImanager* UI = G4UImanager::GetUIpointer();

    // General Menu :
    UI->ApplyCommand("/gui/addMenu file File");
    //        UI->ApplyCommand("/gui/addButton file Continue continue");
    UI->ApplyCommand("/gui/addButton file Exit \"exit\"");

    UI->ApplyCommand("/gui/addMenu run Run");
    UI->ApplyCommand("/gui/addButton run \"beamOn 1\" \"/run/beamOn 1\"");
    UI->ApplyCommand("/gui/addButton run \"beamOn ...\" \"/run/beamOn\"");

    UI->ApplyCommand("/gui/addMenu gun Gun");
    UI->ApplyCommand("/gui/addButton gun \"setDefault\" \"/gun/setDefault\"");
    UI->ApplyCommand("/gui/addButton gun \"Select particle\"  \"/gun/particle\"");
    UI->ApplyCommand("/gui/addButton gun \"Select energy\"  \"/gun/energy\"");
    UI->ApplyCommand("/gui/addButton gun \"Select position\"  \"/gun/position\"");
    UI->ApplyCommand("/gui/addButton gun \"Select direction\"  \"/gun/direction\"");

    // Viewer menu :
    UI->ApplyCommand("/gui/addMenu viewer Viewer");
    UI->ApplyCommand("/gui/addButton viewer \"Set style surface\" \"/vis/viewer/set/style surface\"");
    UI->ApplyCommand("/gui/addButton viewer \"Set style wireframe\" \"/vis/viewer/set/style wire\"");
    UI->ApplyCommand("/gui/addButton viewer \"Refresh viewer\" \"/vis/viewer/refresh\"");
    UI->ApplyCommand("/gui/addButton viewer \"Update viewer (interaction or end-of-file)\" \"/vis/viewer/update\"");
    UI->ApplyCommand("/gui/addButton viewer \"Flush viewer (= refresh + update)\" \"/vis/viewer/flush\"");
    UI->ApplyCommand("/gui/addButton viewer \"Update scene\" \"/vis/scene/notifyHandlers\"");

    // Les commandes icones ne fonctionnent qu avec Qt - pas encore implementees dans autres GUI
    // open/save icons
    UI->ApplyCommand("/gui/addIcon \"Open macro file\" open /control/execute");
    UI->ApplyCommand("/gui/addIcon \"Save viewer state\" save /vis/viewer/save");

    // Cursors style icons
    UI->ApplyCommand("/gui/addIcon Move move");
    UI->ApplyCommand("/gui/addIcon Pick pick");
    UI->ApplyCommand("/gui/addIcon \"Zoom out\" zoom_out");
    UI->ApplyCommand("/gui/addIcon \"Zoom in\" zoom_in");
    UI->ApplyCommand("/gui/addIcon Rotate rotate");

    // Surface Style icons
    UI->ApplyCommand("/gui/addIcon \"Hidden line removal\" hidden_line_removal");
    UI->ApplyCommand("/gui/addIcon \"Hidden line and hidden surface removal\" hidden_line_and_surface_removal");
    UI->ApplyCommand("/gui/addIcon Surfaces solid");
    UI->ApplyCommand("/gui/addIcon Wireframe wireframe");

    // Perspective/Ortho icons
    UI->ApplyCommand("/gui/addIcon Perspective perspective");
    UI->ApplyCommand("/gui/addIcon Orthographic ortho");
  }
#endif
}
