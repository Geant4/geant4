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
// $Id: TargetEventAction.cc,v 1.1 2003-10-08 12:32:11 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "TargetEventAction.hh"

#include "TargetHit.hh"
#include "TargetEventActionMessenger.hh"

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

#include "G4RunManager.hh"
#include "G4Run.hh"


TargetEventAction::TargetEventAction()
:targetCollID(-1),drawFlag("all"),printModulo(1),
 eventMessenger(0),
 ThetaHistP(0), KEHistP(0),
 ThetaHistM(0), KEHistM(0),
 ThetaHistPip30MeV(0),  ThetaHistPip52MeV(0),  ThetaHistPip79MeV(0), 
 ThetaHistPip105MeV(0), ThetaHistPip155MeV(0), ThetaHistPip205MeV(0), 
 ThetaHistPip255MeV(0), ThetaHistPip305MeV(0), ThetaHistPip358MeV(0), 
 ThetaHistPip408MeV(0), ThetaHistPip486MeV(0), ThetaHistPip553MeV(0),
 
 ThetaHistPim30MeV(0),  ThetaHistPim52MeV(0),  ThetaHistPim79MeV(0), 
 ThetaHistPim105MeV(0), ThetaHistPim155MeV(0), ThetaHistPim205MeV(0), 
 ThetaHistPim255MeV(0), ThetaHistPim305MeV(0), ThetaHistPim358MeV(0), 
 ThetaHistPim408MeV(0), ThetaHistPim486MeV(0), ThetaHistPim553MeV(0)
{

  eventMessenger = new TargetEventActionMessenger(this);

// JAIDA start
// 

  //   plotter = 0;
  //   tuple = 0;

  analysisManager = TargetAnalysisManager::getInstance();
  hFactory = analysisManager->getHistogramFactory();

//   ITupleFactory* tFactory = analysisManager->getTupleFactory();
//
//   if (tFactory) {
//      tuple = tFactory->create("Tuple","Tuple","int RunID, double depositEnergy","");
//   }
// 
// JAIDA end

}


TargetEventAction::~TargetEventAction()
{
// JAIDA start
//
   TargetAnalysisManager::dispose();
//
// JAIDA end

  delete eventMessenger;
}


void TargetEventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    HepRandom::showEngineStatus();
  }

  if (targetCollID==-1) {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    targetCollID = SDman->GetCollectionID("TgtCollection");
  }
 
// JAIDA

  if( evtNb == 0 && hFactory ) {  
    G4RunManager* RM = G4RunManager::GetRunManager();
    const G4Run* CurrentRun = RM->GetCurrentRun();
    RunID = CurrentRun -> GetRunID();

    G4String HistName;
    HistName = "Run #";
    HistName = HistName + itoa ( RunID );
    ThetaHistP = hFactory->create1D("Pi+ cos(theta)",200,-1.0, 1.0);
    KEHistP = hFactory->create1D("Pi+ KE(MeV)",100,0.,1000.0);
    ThetaHistPip30MeV = hFactory->create1D("Pi+ 30 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip52MeV = hFactory->create1D("Pi+ 52 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip79MeV = hFactory->create1D("Pi+ 79 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip105MeV = hFactory->create1D("Pi+ 105 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip155MeV = hFactory->create1D("Pi+ 155 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip205MeV = hFactory->create1D("Pi+ 205 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip255MeV = hFactory->create1D("Pi+ 255 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip305MeV = hFactory->create1D("Pi+ 305 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip358MeV = hFactory->create1D("Pi+ 358 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip408MeV = hFactory->create1D("Pi+ 408 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip486MeV = hFactory->create1D("Pi+ 486 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPip553MeV = hFactory->create1D("Pi+ 553 MeV cos(theta)",
                                                            100,-1.0,1.0);

    ThetaHistM = hFactory->create1D("Pi- cos(theta)",200,-1.0, 1.0);
    KEHistM = hFactory->create1D("Pi- KE(MeV)",100,0.,1000.0);
    ThetaHistPim30MeV = hFactory->create1D("Pi- 30 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim52MeV = hFactory->create1D("Pi- 52 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim79MeV = hFactory->create1D("Pi- 79 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim105MeV = hFactory->create1D("Pi- 105 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim155MeV = hFactory->create1D("Pi- 155 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim205MeV = hFactory->create1D("Pi- 205 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim255MeV = hFactory->create1D("Pi- 255 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim305MeV = hFactory->create1D("Pi- 305 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim358MeV = hFactory->create1D("Pi- 358 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim408MeV = hFactory->create1D("Pi- 408 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim486MeV = hFactory->create1D("Pi- 486 MeV cos(theta)",
                                                            100,-1.0,1.0);
    ThetaHistPim553MeV = hFactory->create1D("Pi- 553 MeV cos(theta)",
                                                            100,-1.0,1.0);

//      plotter = analysisManager->createPlotter();
//      plotter->createRegions(1,1);
//      plotter->plot(*TargHist);
//      plotter->show();
  }
// JAIDA end
}


void TargetEventAction::EndOfEventAction(const G4Event* evt)
{
  //  G4int evtNb = evt->GetEventID();

  // Extract info from hits

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  TargetHitsCollection* THC = 0;

  G4int n_hit = 0;
  G4double cosTheta = 0;
  G4double kineticEnergy = 0;
  G4double charge = 0;

  if (HCE) THC = (TargetHitsCollection*)(HCE->GetHC(targetCollID));

  if (THC) {
     n_hit = THC->entries();
     for (G4int i=0;i<n_hit;i++) {
       cosTheta = (*THC)[i]->GetTheta(); 
       kineticEnergy = (*THC)[i]->GetEkin();
       charge = (*THC)[i]->GetCharge();
     }
  }

// JAIDA start
//
  if(charge > 0.0) {
    // pi+ histograms
    //
    ThetaHistP->fill( cosTheta );
    KEHistP->fill( kineticEnergy );

    if(kineticEnergy > 27.75 && kineticEnergy < 32.30) {
      ThetaHistPip30MeV->fill(cosTheta);
    } else if(kineticEnergy > 47.58 && kineticEnergy < 56.48) { 
      ThetaHistPip52MeV->fill(cosTheta);
    } else if(kineticEnergy > 73.27 && kineticEnergy < 84.81) { 
      ThetaHistPip79MeV->fill(cosTheta);
    } else if(kineticEnergy > 98.11 && kineticEnergy < 112.07) { 
      ThetaHistPip105MeV->fill(cosTheta);
    } else if(kineticEnergy > 146.21 && kineticEnergy < 163.82) { 
      ThetaHistPip155MeV->fill(cosTheta);
    } else if(kineticEnergy > 195.01 && kineticEnergy < 215.12) { 
      ThetaHistPip205MeV->fill(cosTheta);
    } else if(kineticEnergy > 243.82 && kineticEnergy < 266.26) { 
      ThetaHistPip255MeV->fill(cosTheta);
    } else if(kineticEnergy > 292.66 && kineticEnergy < 317.35) { 
      ThetaHistPip305MeV->fill(cosTheta);
    } else if(kineticEnergy > 345.04 && kineticEnergy < 370.96) { 
      ThetaHistPip358MeV->fill(cosTheta);
    } else if(kineticEnergy > 393.99 && kineticEnergy < 422.03) { 
      ThetaHistPip408MeV->fill(cosTheta);
    } else if(kineticEnergy > 471.36 && kineticEnergy < 500.61) { 
      ThetaHistPip486MeV->fill(cosTheta);
    } else if(kineticEnergy > 537.84 && kineticEnergy < 568.20) { 
      ThetaHistPip553MeV->fill(cosTheta);
    }
  } else if (charge < 0.0) {
    // pi- histograms
    //
    ThetaHistM->fill( cosTheta );
    KEHistM->fill( kineticEnergy );

    if(kineticEnergy > 27.75 && kineticEnergy < 32.30) { 
      ThetaHistPim30MeV->fill(cosTheta);
    } else if(kineticEnergy > 47.58 && kineticEnergy < 56.48) { 
      ThetaHistPim52MeV->fill(cosTheta);
    } else if(kineticEnergy > 73.27 && kineticEnergy < 84.81) { 
      ThetaHistPim79MeV->fill(cosTheta);
    } else if(kineticEnergy > 98.11 && kineticEnergy < 112.07) { 
      ThetaHistPim105MeV->fill(cosTheta);
    } else if(kineticEnergy > 146.21 && kineticEnergy < 163.82) { 
      ThetaHistPim155MeV->fill(cosTheta);
    } else if(kineticEnergy > 195.01 && kineticEnergy < 215.12) { 
      ThetaHistPim205MeV->fill(cosTheta);
    } else if(kineticEnergy > 243.82 && kineticEnergy < 266.26) { 
      ThetaHistPim255MeV->fill(cosTheta);
    } else if(kineticEnergy > 292.66 && kineticEnergy < 317.35) { 
      ThetaHistPim305MeV->fill(cosTheta);
    } else if(kineticEnergy > 345.04 && kineticEnergy < 370.96) { 
      ThetaHistPim358MeV->fill(cosTheta);
    } else if(kineticEnergy > 393.99 && kineticEnergy < 422.03) { 
      ThetaHistPim408MeV->fill(cosTheta);
    } else if(kineticEnergy > 471.36 && kineticEnergy < 500.61) { 
      ThetaHistPim486MeV->fill(cosTheta);
    } else if(kineticEnergy > 537.84 && kineticEnergy < 568.20) { 
      ThetaHistPim553MeV->fill(cosTheta);
    }
  }

  /* 
      
//   if (tuple)
//   {
//      tuple->fill( 0 , RunID );
//      tuple->fill( 1 , totEAbs + totEGap );
//      tuple->addRow();

//
// JAIDA end


  // extract the trajectories and draw them
  
  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                                ((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") trj->DrawTrajectory(50);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(50);
          else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
                                  trj->DrawTrajectory(50);
        }
  }
  */
} 









