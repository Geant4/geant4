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
//
// G4VRMLView.cc
// Yasuhide Sawada & Satoshi Tanaka

#ifndef WIN32

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include "G4VisManager.hh"

#include "G4Scene.hh"
#include "G4VRML1Viewer.hh"
#include "G4VRML1SceneHandler.hh"
#include "G4VRML1.hh"
#include "G4ios.hh"

G4VRML1Viewer::G4VRML1Viewer(G4VRML1SceneHandler& scene, const G4String& name) :
   G4VViewer(scene, scene.IncrementViewCount(), name), fSceneHandler(scene)
{}

G4VRML1Viewer::~G4VRML1Viewer()
{}

void G4VRML1Viewer::SetView()
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
        G4cout << "***** G4VRML1Viewer::SetView(): No effects" << G4endl;
#endif
}

void G4VRML1Viewer::DrawView()
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4VRML1Viewer::DrawView()" << G4endl;
#endif

	fSceneHandler.VRMLBeginModeling();

	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML1Viewer::ClearView(void)
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
        G4cout << "***** G4VRML1Viewer::ClearView(): No effects" << G4endl;
#endif
}

void G4VRML1Viewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
        G4cout << "***** G4VRML1Viewer::ShowView()" << G4endl;
#endif
	fSceneHandler.VRMLEndModeling();
}

void G4VRML1Viewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
        G4cout << "***** G4VRML1Viewer::FinishView(): No effects" << G4endl;
#endif
}

#endif
#endif //WIN32
