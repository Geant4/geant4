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

#include "G4TheMTRayTracer.hh"
#include "G4SystemOfUnits.hh"
#include "G4RTMessenger.hh"
#include "G4VFigureFileMaker.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTRun.hh"
#include "G4RTRunAction.hh"
#include "G4RTWorkerInitialization.hh"
#include "G4VRTScanner.hh"

#include "G4MTRunManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4VVisManager.hh"
#include "G4RunManagerFactory.hh"

#define G4warn G4cout

G4TheMTRayTracer* G4TheMTRayTracer::theInstance = nullptr;

G4TheMTRayTracer::G4TheMTRayTracer(G4VFigureFileMaker* figMaker,
			       G4VRTScanner* scanner)
: G4TheRayTracer(figMaker,scanner)
{
  if(!theInstance) 
  { theInstance = this; }
  else
  { G4Exception("G4TheMTRayTracer::G4TheMTRayTracer","VisRayTracer00100",
         FatalException,"G4TheMTRayTracer has to be a singleton.");}
  theUserWorkerInitialization = 0;
  theRTWorkerInitialization = 0;
  theUserRunAction = 0;
  theRTRunAction = 0;
}

G4TheMTRayTracer* G4TheMTRayTracer::Instance()
{
  if (theInstance) return theInstance;
  else return new G4TheMTRayTracer;
}

G4TheMTRayTracer* G4TheMTRayTracer::Instance
(G4VFigureFileMaker* figMaker,G4VRTScanner* scanner)
{
  if (theInstance) {
    theFigMaker=figMaker;
    theScanner=scanner;
    return theInstance;
  }
  else return new G4TheMTRayTracer(figMaker,scanner);
}

G4TheMTRayTracer::~G4TheMTRayTracer()
{
  if(theRTWorkerInitialization)
  {
    delete theRTWorkerInitialization;
    theRTWorkerInitialization = 0;
  }
  if(theRTRunAction)
  {
    delete theRTRunAction;
    theRTRunAction = 0;
  }
}

void G4TheMTRayTracer::Trace(const G4String& fileName)
{
  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  G4ApplicationState currentState = theStateMan->GetCurrentState();
  if(currentState!=G4State_Idle)
  {
    G4warn << "Illegal application state <" << theStateMan->GetStateString(currentState)
           << "> - Trace() ignored. " << G4endl;
    return;
  }

  if(!theFigMaker)
  {
    G4warn << "Figure file maker class is not specified - Trace() ignored." << G4endl;
    return;
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int storeTrajectory = UI->GetCurrentIntValue("/tracking/storeTrajectory");
  UI->ApplyCommand("/tracking/storeTrajectory 1");

  G4ThreeVector tmpVec = targetPosition - eyePosition;
  eyeDirection = tmpVec.unit();
  G4int nPixel = nColumn*nRow;
  colorR = new unsigned char[nPixel];
  colorG = new unsigned char[nPixel];
  colorB = new unsigned char[nPixel];
  unsigned char defR = (unsigned char)(G4int(255*backgroundColour.GetRed()));
  unsigned char defG = (unsigned char)(G4int(255*backgroundColour.GetGreen()));
  unsigned char defB = (unsigned char)(G4int(255*backgroundColour.GetBlue()));
  for(G4int ii=0;ii<nPixel;++ii)
  {
    colorR[ii] = defR;
    colorG[ii] = defG;
    colorB[ii] = defB;
  }

  G4bool succeeded = CreateBitMap();
  if(succeeded)
  { CreateFigureFile(fileName); }
  else
  { G4warn << "Could not create figure file" << G4endl;
    G4warn << "You might set the eye position outside of the world volume" << G4endl; }

  G4String str = "/tracking/storeTrajectory " + G4UIcommand::ConvertToString(storeTrajectory);
  UI->ApplyCommand(str);

  delete [] colorR;
  delete [] colorG;
  delete [] colorB;
}

void G4TheMTRayTracer::StoreUserActions()
{
  G4MTRunManager* mrm         = G4RunManagerFactory::GetMTMasterRunManager();
  theUserWorkerInitialization = mrm->GetUserWorkerInitialization();
  theUserRunAction = mrm->GetUserRunAction();

  if(!theRTWorkerInitialization) theRTWorkerInitialization = new G4RTWorkerInitialization();
  if(!theRTRunAction) theRTRunAction = new G4RTRunAction();

  mrm->SetUserInitialization(theRTWorkerInitialization);
  mrm->SetUserAction(theRTRunAction);
}

void G4TheMTRayTracer::RestoreUserActions()
{
  G4MTRunManager* mrm = G4RunManagerFactory::GetMTMasterRunManager();
  mrm->SetUserInitialization(
    const_cast<G4UserWorkerInitialization*>(theUserWorkerInitialization));
  mrm->SetUserAction(const_cast<G4UserRunAction*>(theUserRunAction));
}

G4bool G4TheMTRayTracer::CreateBitMap()
{
  G4VVisManager* visMan = G4VVisManager::GetConcreteInstance();
  visMan->IgnoreStateChanges(true);
  StoreUserActions();

  G4MTRunManager* mrm = G4RunManagerFactory::GetMTMasterRunManager();

  // Keep, then switch off any printing requests
  auto runVerbosity      = mrm->GetVerboseLevel();
  auto runPrintProgress  = mrm->GetPrintProgress();
  G4UImanager::GetUIpointer()->ApplyCommand("/run/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/run/printProgress 0");

  // Event loop
  G4int nEvent = nRow*nColumn;
////  mrm->BeamOn(nEvent);
////  Temporary work-around until direct invokation of G4RunManager::BeamOn() works.
  G4String str = "/run/beamOn " + G4UIcommand::ConvertToString(nEvent);
  G4UImanager::GetUIpointer()->ApplyCommand(str);

  // Restore printing requests
  str = "/run/verbose " + G4UIcommand::ConvertToString(runVerbosity);
  G4UImanager::GetUIpointer()->ApplyCommand(str);
  str = "/run/printProgress " + G4UIcommand::ConvertToString(runPrintProgress);
  G4UImanager::GetUIpointer()->ApplyCommand(str);

  RestoreUserActions();
  visMan->IgnoreStateChanges(false);

  const G4RTRun* theRun = static_cast<const G4RTRun*>(mrm->GetCurrentRun());
  if(!theRun) return false;

  G4THitsMap<G4Colour>* colMap = theRun->GetMap(); 
  auto itr = colMap->GetMap()->cbegin();
  for(;itr!=colMap->GetMap()->cend();++itr)
  {
    G4int key = itr->first;
    G4Colour* col = itr->second;
    colorR[key] = (unsigned char)(G4int(255*col->GetRed()));
    colorG[key] = (unsigned char)(G4int(255*col->GetGreen()));
    colorB[key] = (unsigned char)(G4int(255*col->GetBlue()));
  }  

  theScanner->Initialize(nRow,nColumn);
  G4int iRow, iColumn;
  while (theScanner->Coords(iRow,iColumn))
  {
    G4int iCoord = iRow * nColumn + iColumn;
    theScanner->Draw(colorR[iCoord],colorG[iCoord],colorB[iCoord]);
  }

  return true;
}

