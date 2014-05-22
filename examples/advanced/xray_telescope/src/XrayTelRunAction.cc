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

#include "XrayTelRunAction.hh"
#include "XrayTelAnalysis.hh"

#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

XrayTelRunAction::XrayTelRunAction()
{ }


XrayTelRunAction::~XrayTelRunAction()
{ }


void XrayTelRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int runN = aRun->GetRunID();
  if (IsMaster())
    G4cout << "### Run : " << runN << " (master)" << G4endl;
  else
    G4cout << "### Run : " << runN << " (worker)" << G4endl;

  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
  // Book histograms and ntuples
  XrayTelAnalysis* analysis = XrayTelAnalysis::getInstance();
  analysis->book(IsMaster());
}


void XrayTelRunAction::EndOfRunAction(const G4Run* )
{
  XrayTelAnalysis* analysis = XrayTelAnalysis::getInstance();
  analysis->finish(IsMaster());

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
}





















