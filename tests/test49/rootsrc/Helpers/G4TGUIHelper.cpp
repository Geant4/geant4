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
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TGUIHelper.h"

G4TGUIHelper *gGUIHelper = new G4TGUIHelper();


ClassImp(G4TGUIHelper)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
void G4TGUIHelper::ShowMenu()
{
	G4TSimHelper::LoadLibraries();
	G4TToolMenuGUI* Window = new G4TToolMenuGUI(gClient->GetRoot(), 0, 0);
	Window->MapWindow();
}

//______________________________________________________________________________
void G4TGUIHelper::ShowSimulationGUI()
{
	G4TSimHelper::LoadLibraries();
	G4TSimulationGUI* Window = new G4TSimulationGUI(gClient->GetRoot(), 0, 0);
	Window->MapWindow();
}


//______________________________________________________________________________
void G4TGUIHelper::ShowAnalysisGUI()
{
	G4TSimHelper::LoadLibraries();
	G4TAnalysisGUI* Window = new G4TAnalysisGUI(gClient->GetRoot(), 0, 0);
	Window->MapWindow();
}

//______________________________________________________________________________
void G4TGUIHelper::ShowPublicationGUI()
{
	G4TSimHelper::LoadLibraries();
	G4TPublicationGUI* Window = new G4TPublicationGUI(gClient->GetRoot(), 0, 0);
	Window->MapWindow();
}


