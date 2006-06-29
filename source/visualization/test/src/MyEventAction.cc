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
// $Id: MyEventAction.cc,v 1.14 2006-06-29 21:34:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyEventAction.hh"

//#define DRAWTRAJHIT

#ifdef DRAWTRAJHIT

#include "MyTrackerHit.hh"
#include "MyCalorimeterHit.hh"

#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Scale.hh"
#include "G4Text.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include <sstream>

MyEventAction::MyEventAction()
{;}

MyEventAction::~MyEventAction()
{;}

void MyEventAction::BeginOfEventAction(const G4Event*)
{;}

void MyEventAction::EndOfEventAction(const G4Event* anEvent)
{
  static int coutCount = 0;
  if (coutCount < 10) {
    coutCount++;
    G4cout << "MyEventAction::EndOfEventActionAction called." << G4endl;
  }

  const G4Event* evt = anEvent;

#ifdef DRAWTRAJHIT

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int trackerCollID = SDman->GetCollectionID(colNam="TrackerCollection");
  G4int calorimeterCollID = SDman->GetCollectionID(colNam="CalCollection");

#endif

  if (coutCount < 10) {
    G4cout << ">>> Run "
	<< G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID()
	<< " Event " << anEvent->GetEventID() << G4endl;
  }

#ifdef DRAWTRAJHIT

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
  { n_trajectories = trajectoryContainer->entries(); }
  if (coutCount < 10) {
    G4cout << "    " << n_trajectories 
	   << " trajectories stored in this event." << G4endl;
  }

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  G4int n_hitCollection = 0;
  if(HCE)
  { n_hitCollection = HCE->GetCapacity(); }
  if (coutCount < 10) {
    G4cout << "    " << n_hitCollection
	   << " hitsCollections stored in this event." << G4endl;
  }

#endif

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {

#ifdef DRAWTRAJHIT

    for(G4int i=0; i<n_trajectories; i++)
    { (*(evt->GetTrajectoryContainer()))[i]->DrawTrajectory(50); }

    MyTrackerHitsCollection* THC 
      = (MyTrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    if(THC) THC->DrawAllHits();
    MyCalorimeterHitsCollection* CHC
      = (MyCalorimeterHitsCollection*)(HCE->GetHC(calorimeterCollID));
    if(CHC) CHC->DrawAllHits();

#endif

    G4Scale scale(1. * m, "Test Scale");
    G4VisAttributes va (G4Colour(1.,0.,0.));
    scale.SetVisAttributes(va);
    pVVisManager->Draw(scale);

    std::ostringstream oss;
    oss << "Run "
	<< G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID()
	<< " Event " << anEvent->GetEventID();
    //G4Text text(oss.str(), G4Point3D(400.*cm, 400.*cm, -400.*cm));
    G4Text text(oss.str(), G4Point3D(0., -0.9, 0.));
    text.SetScreenSize(18);
    G4VisAttributes textAtts(G4Colour(0.,1.,1));
    text.SetVisAttributes(textAtts);
    pVVisManager->Draw2D(text);

    G4Box transientBox("transientBox",100*cm,100*cm,100*cm);
    G4VisAttributes transientBoxAtts(G4Colour(1.,0.,1));
    transientBoxAtts.SetForceWireframe(true);
    pVVisManager->Draw(transientBox, transientBoxAtts,
		       G4Translate3D(500.*cm, 500.*cm, -500.*cm));

    G4Tubs transientTube("transientTube",0.,100*cm,100*cm,0.,360.*deg);
    G4VisAttributes transientTubeAtts(G4Colour(1.,1.,0));
    transientTubeAtts.SetForceWireframe(true);
    pVVisManager->Draw(transientTube, transientTubeAtts,
		       G4Translate3D(500.*cm, 300.*cm, -500.*cm));

  }
}
