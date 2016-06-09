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
// Rich advanced example for Geant4
// RichTbAnalysisManager.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "G4Timer.hh"
#include "globals.hh"
#include "fstream"
#include <ctype.h> 
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4SteppingManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "RichTbHit.hh"
#include "RichTbRunConfig.hh"

#include "RichTbAnalysisMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "RichTbAnalysisManager.hh"
#include "AIDA/AIDA.h"



RichTbAnalysisManager* RichTbAnalysisManager::instance = 0;

RichTbAnalysisManager::RichTbAnalysisManager()
  :outputFileName("rich.his"), analysisFactory(0),tree(0),histogramFactory(0)
{

  iTimer = new G4Timer;
  NumPhotBeforeAerogel=0;
  NumPhotBeforeMirror=0;
  NumPhotAfterMirror=0;
  NumPhotBeforeFilter=0;
  NumPhotAfterFilter=0;
  NumPhotAtHpd1Input=0;
  NumPhotAtHpd2Input=0;
  NumPhotAtHpd3Input=0;
  NumPhotAtHpd4Input=0;
  NumHitTotInHpd1=0;
  NumHitTotInHpd2=0;
  NumHitTotInHpd3=0;
  NumHitTotInHpd4=0;
  NumHitInSi=0;


  analisysMessenger = new RichTbAnalysisMessenger(this);

// Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {

    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {

      tree = treeFactory->create(outputFileName,"hbook",false,true);

      delete treeFactory; // Will not delete the ITree.
      histogramFactory = analysisFactory->createHistogramFactory(*tree);  

    }
  }

}

RichTbAnalysisManager::~RichTbAnalysisManager(){
  
  delete histogramFactory;
  histogramFactory = 0;

  delete analysisFactory;
  analysisFactory = 0;

  delete instance;
}

RichTbAnalysisManager* RichTbAnalysisManager::getInstance()
{
  if(instance==0) {instance = new RichTbAnalysisManager;}
  return instance;
}

void RichTbAnalysisManager::book()
{

  fhistoNrPhotG = histogramFactory->createHistogram1D("1","Number of Photon Hits per Track",120,0.,20.);
  fhistoNBeforeMirror = histogramFactory->createHistogram1D("2","Number of Photons before Mirror",100,0.,750.); 
  fhistoWBeforeMirror = histogramFactory->createHistogram1D("3","WaveLength of Photons before Mirror",100,200.,800.);
  fhistoWAfterMirror = histogramFactory->createHistogram1D("4","WaveLength of Photons after Mirror",100,200.,800.);
  fhistoCkvProdSmall = histogramFactory->createHistogram1D("9","Cherekov Angle at roduction all angles",500,0.,0.5);
  fhistoEmisZ = histogramFactory->createHistogram1D("10","Z of the Photon Emission Point",120,0.,600.);
  fhistoCkvRadius = histogramFactory->createHistogram1D("5","Cherenkov Ring Radius",300,80.,220.);

}

void RichTbAnalysisManager::finish()
{
  if(tree){

    tree->commit(); //Write histos and ntuples in file
    tree->close();
  }

}


void RichTbAnalysisManager::BeginOfEventAnalysis(const G4Event*){

  //RichCollId is already defined in LHCbRichSimEventAction.cc
  // Hence its extraction is not repeated here.
  iTimer->Start();
  NumPhotBeforeAerogel=0;
  NumPhotBeforeMirror=0;
  NumPhotAfterMirror=0;
  NumPhotBeforeFilter=0;
  NumPhotAfterFilter=0;
  NumPhotAtHpd1Input=0;
  NumPhotAtHpd2Input=0;
  NumPhotAtHpd3Input=0;
  NumPhotAtHpd4Input=0;
  NumHitTotInHpd1=0;
  NumHitTotInHpd2=0;
  NumHitTotInHpd3=0;
  NumHitTotInHpd4=0;
  NumHitInSi=0;
  
}
void RichTbAnalysisManager::EndOfEventAnalysis(const G4Event* evt){



  iTimer->Stop();
  // G4double TimeforThisEvent= iTimer->GetRealElapsed();

  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  G4String colNam;
  G4int RichTbCollID = SDman->GetCollectionID(colNam="RichTbHitsCollection");

  G4int RichTbc =  RichTbCollID;
  G4int evtNb = evt->GetEventID();

  if(RichTbc<0) return;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  RichTbHitsCollection* RHC = 
    HCE ? (RichTbHitsCollection*)(HCE->GetHC(RichTbc))  : 0;
  if(RHC)
   {

    int n_hit = RHC->entries();
    G4cout << " " << n_hit
      << " hits are stored in RichTbHitsCollection in event "<<evtNb<< G4endl;
    //fill the final number of hits in this event into a histogram

    fhistoNrPhotG->fill(n_hit*1.0);

   //To print and plot the number of hits per each HPD.

     for (G4int iha=0; iha<n_hit; iha++ ) {

      RichTbHit* aHit = (*RHC)[iha];
      G4int CurHpdnum = aHit -> GetCurHpdNum();
      if(CurHpdnum == 0 ) {
	NumHitTotInHpd1++;
      
      }else if (CurHpdnum == 1 ) {
	NumHitTotInHpd2++;

      }else if (CurHpdnum == 2 ) {
	NumHitTotInHpd3++;
      
      }else if (CurHpdnum == 3 ) {
	NumHitTotInHpd4++;

      }

     }

    //fill the number of photons before mirror into a histogram

     fhistoNBeforeMirror->fill(NumPhotBeforeMirror);

     G4cout<<"NumPhot before mirror "<<NumPhotBeforeMirror<<G4endl;
     G4cout<<"NumPhot after mirror "<<NumPhotAfterMirror<<G4endl;
     G4cout<<"NumPhot incident on Hpd 1 2 3 4 = "<< NumPhotAtHpd1Input
	   <<"  "<<  NumPhotAtHpd2Input<<"  "<< NumPhotAtHpd3Input
	   <<"  "<< NumPhotAtHpd4Input<<G4endl;
     
     G4cout<<"Number of Si hits in Hpd 1 2 3 4 are   "<<NumHitTotInHpd1
           <<"  "<<NumHitTotInHpd2<<"   "<<NumHitTotInHpd3
           <<"  "<<NumHitTotInHpd4<<G4endl;
     
     
   }
  
}
//void RichTbAnalysisManager::StepAnalysis(
//       const G4SteppingManager* aSteppingManager ){

//}
void  RichTbAnalysisManager::setMeanHXCoord(G4double cx ) {

  
  G4double hcurx = NumHXCoord+1.0;
  if(hcurx != 0.0 ) {
  G4double Mhx =  MeanHXCoord * (NumHXCoord / hcurx) + cx/hcurx ;
  MeanHXCoord = Mhx;
  NumHXCoord  = hcurx;

  }else { G4cout<< "Zero hcurx in setMeanHXCoord "<<G4endl; }
}
void  RichTbAnalysisManager::setMeanHYCoord(G4double cy ) {

 G4double hcury = NumHYCoord+1.0;
  if(hcury != 0.0 ) {
  G4double Mhy =  MeanHYCoord * (NumHYCoord / hcury) + cy/hcury ;
  MeanHYCoord = Mhy;
  NumHYCoord  = hcury;
   }else { G4cout<< "Zero hcury in setMeanHYCoord "<<G4endl; } 
}


void RichTbAnalysisManager::SetOutputFileName(G4String newName)
{
  outputFileName = newName;
}








#endif

