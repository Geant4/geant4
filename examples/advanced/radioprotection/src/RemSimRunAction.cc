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
// ********************************************************************

//
// --------------------------------------------------------------
//                 GEANT 4 - RemSimtherapy example
// --------------------------------------------------------------
//
// Code developed by:
//  S.Guatelli
//
//
//    *******************************
//    *                             *
//    *    RemSimRunAction.cc       *
//    *                             *
//    *******************************
//
// $Id: RemSimRunAction.cc,v 1.3 2004-03-12 10:55:55 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimRunAction.hh"
#include "RemSimDetectorConstruction.hh"

#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "RemSimRunAction.hh"
#include <fstream>
#include <strstream>

RemSimRunAction::RemSimRunAction()
{
}

RemSimRunAction::~RemSimRunAction()
{   
 }
void RemSimRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
 runID = aRun->GetRunID();
}

void RemSimRunAction::EndOfRunAction(const G4Run* aRun)
{  
G4double numberEvents = aRun -> GetNumberOfEvent();
 G4cout<<"Number of events: "<<numberEvents<<G4endl;
}
