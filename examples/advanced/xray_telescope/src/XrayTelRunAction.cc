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
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelRunAction.cc                             *
// * -------                                                            *
// *                                                                    *
// * Version:           0.5                                             *
// * Date:              16/10/01                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 07.11.2001 M.G. Pia
// - Modified the analysis management
// - Small design iteration
//
// 16.10.2001 R. Nartallo
// - Updated "/vis" commands to new versions
// - Clean up code to avoid 'pedantic' and 'ANSI' compiler warnings 
//
// 30.11.2000 R. Nartallo
// - Add pre-processor directives to compile without analysis option
//
// 16.11.2000 A. Pfeiffer
// - Implementation of analysis manager call
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope Physics list
// - Based on Chandra and XMM models
// 
//
// **********************************************************************

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "XrayTelRunAction.hh"
#include "XrayTelAnalysis.hh"
#include "SystemOfUnits.h"


XrayTelRunAction::XrayTelRunAction()
  :nEnteringTracks(0), totEnteringEnergy(0.)
{ }


XrayTelRunAction::~XrayTelRunAction()
{ }


void XrayTelRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int runN = aRun->GetRunID();
  if ( runN % 1000 == 0 ) 
    G4cout << "### Run : " << runN << G4endl;

  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 

  nEnteringTracks = 0;
  totEnteringEnergy = 0.;

  // Book histograms and ntuples
  XrayTelAnalysis* analysis = XrayTelAnalysis::getInstance();
  analysis->book();
}


void XrayTelRunAction::EndOfRunAction(const G4Run* )
{
  XrayTelAnalysis* analysis = XrayTelAnalysis::getInstance();
  analysis->finish();

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  G4cout << "End of Run summary" << G4endl << G4endl;

  G4cout << "Total Entering Detector : " << nEnteringTracks  << G4endl;
  G4cout << "Total Entering Detector Energy : " 
	 << totEnteringEnergy/MeV  
	 << " MeV"
	 << G4endl;
}


void XrayTelRunAction::Update(G4double energy)
{
  nEnteringTracks++;
  totEnteringEnergy += energy;
}
























