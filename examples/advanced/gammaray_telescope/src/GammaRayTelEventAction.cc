// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelEventAction.cc,v 1.2 2000-11-15 20:27:41 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelEventAction  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelEventAction.hh"
#include "GammaRayTelPayloadHit.hh"
#include "g4rw/tvordvec.h"

#ifdef G4HIS_USE_AIDA
#include "GammaRayTelHistogram.hh"
#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

// This file is a global variabe in which we store energy deposition per hit
// and other relevant information
extern G4std::ofstream outFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4HIS_USE_AIDA
GammaRayTelEventAction::GammaRayTelEventAction(GammaRayTelHistogram* hMgr)
  :drawFlag("all"), printModulo(10000),trackerCollID(-1),histoManager(hMgr)
{
}
#else
GammaRayTelEventAction::GammaRayTelEventAction()
  :drawFlag("all"), printModulo(10000),trackerCollID(-1)
{
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::~GammaRayTelEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0)
    { 
      G4cout << "\n---> Event: " << evtNb << G4endl;
    }
  
  if (trackerCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      trackerCollID = SDman->GetCollectionID("PayloadCollection");
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  
#ifdef G4HIS_USE_AIDA
  IHistogram1D * Edeposited;
  vector<IHistogram1D *> * hlist = histoManager->getH1DList();
  Edeposited = (*hlist)[0]; // 1D-histo # 0 is kinetic engergy
  IHistogram2D * posXZHist;
  // 2D-histo # 0 is posXZ
  posXZHist = (*(histoManager->getH2DList()))[0];
  IHistogram2D * posYZHist;
  // 2D-histo # 1 is posYZ
  posYZHist = (*(histoManager->getH2DList()))[1];
#endif
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  GammaRayTelPayloadHitsCollection* CHC = NULL;
  if (HCE)
    CHC = (GammaRayTelPayloadHitsCollection*)(HCE->GetHC(trackerCollID));
  
  if (CHC)
    {
      int n_hit = CHC->entries();
      G4cout << "Number of hits in this event =  " << n_hit << G4endl;
      G4double ESil=0;
      G4int NStrip, NPlane, IsX;
      for (int i=0;i<n_hit;i++)
	{
	  // Here I put the energy deposition per hit in the file for 
	  // later analysis 
	  ESil = (*CHC)[i]->GetEdepSil();
	  NStrip = (*CHC)[i]->GetNStrip();
	  NPlane = (*CHC)[i]->GetNSilPlane();
	  IsX = (*CHC)[i]->GetPlaneType();
	  outFile << G4std::setw(7) << ESil/keV << " " << NStrip << 
	    " " << NPlane << " " << IsX << " " <<
	    (*CHC)[i]->GetPos().x() <<" "<<
	    (*CHC)[i]->GetPos().y() <<" "<<
	    (*CHC)[i]->GetPos().z() <<" "<<
	    G4endl;
	  
#ifdef G4HIS_USE_AIDA
	  Edeposited->fill(ESil);
	  if(IsX) 
	    posXZHist->fill((*CHC)[i]->GetPos().x(),
			    (*CHC)[i]->GetPos().z());
	  else
	    posYZHist->fill((*CHC)[i]->GetPos().y(),
			    (*CHC)[i]->GetPos().z());
#endif
	}
      
#ifdef G4HIS_USE_AIDA
      histoManager->plot1(posXZHist);
      histoManager->plot2(posYZHist);
      histoManager->getPlotter1()->refresh();
      histoManager->getPlotter2()->refresh();
      histoManager->getPlotter1()->psPrint("plotXZ.ps");
      histoManager->getPlotter2()->psPrint("plotYZ.ps");
#endif

    }
  
  
  
  if(G4VVisManager::GetConcreteInstance())
    {
      for(G4int i=0; i<n_trajectories; i++) 
	{ G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
	if (drawFlag == "all") trj->DrawTrajectory(50);
	else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	  trj->DrawTrajectory(50); 
	}
    }
}











