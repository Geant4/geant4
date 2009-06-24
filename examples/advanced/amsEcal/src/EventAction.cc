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
// $Id: EventAction.cc,v 1.5 2009-06-24 21:04:00 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventActionMessenger.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det, RunAction* run,
                         PrimaryGeneratorAction* prim, HistoManager* hist)
:detector(det), runAct(run), primary(prim), histoManager(hist)
{ 
  trigger = false;
  Eseuil  = 10*keV;
    
  drawFlag = "none";
  printModulo = 1000;
  eventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{   
  G4int evtNb = evt->GetEventID();

  //survey printing
  if (evtNb%printModulo == 0)
    G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
    
  //initialize Energy per event
  //
  G4int nbOfPixels = detector->GetSizeVectorPixels();
  G4int size = totalEnergy.size();
  if (size < nbOfPixels) {
    visibleEnergy.resize(nbOfPixels);
      totalEnergy.resize(nbOfPixels);
  }
  
  for (G4int k=0; k<nbOfPixels; k++) {
    visibleEnergy[k] = totalEnergy[k] = 0.0;
  }   
  nbRadLen = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int nbOfLayers  = detector->GetNxPixelsTot();
  G4int nyPixels    = detector->GetNyPixels();
  G4int nyPixelsMax = detector->GetNyPixelsMax();
      
  // code for trigger conditions :
  // 1 and only 1 pixel fired per layer
  //
  if (trigger) {
    for (G4int ix=0; ix<nbOfLayers; ix++) {
      //count number of pixels fired
      G4int count = 0;  
      for (G4int iy=0; iy<nyPixels; iy++) {
        G4int k = ix*nyPixelsMax + iy;
        if (visibleEnergy[k] > Eseuil) count++;	      
      }
      //if event killed --> skip EndOfEventAction          
      if (count > 1) return;
    }  
  }
  
  //pass informations to RunAction and HistoManager
  //
  G4double calorEvis = 0.;
  G4double calorEtot = 0.;  
  for (G4int ix=0; ix<nbOfLayers; ix++) {
    //sum energy per readout layer  
    G4double layerEvis = 0.;
    G4double layerEtot = 0.;  
    for (G4int iy=0; iy<nyPixels; iy++) {
      G4int k = ix*nyPixelsMax + iy;
      runAct->fillPerEvent_1(k,visibleEnergy[k],totalEnergy[k]);      
      layerEvis += visibleEnergy[k];
      layerEtot += totalEnergy[k];
      calorEvis += visibleEnergy[k];
      calorEtot += totalEnergy[k];		      
    }      
    runAct->fillPerEvent_2(ix,layerEvis,layerEtot);
    if (layerEvis > 0.) histoManager->FillNtuple(1, ix, layerEvis);
    if (layerEtot > 0.) histoManager->FillNtuple(1, nbOfLayers+ix, layerEtot);
  }
  
  histoManager->AddRowNtuple(1);
  
  if (calorEvis > 0.) histoManager->FillHisto(1,calorEvis);
  if (calorEtot > 0.) histoManager->FillHisto(2,calorEtot);
  
  G4double Ebeam = primary->GetParticleGun()->GetParticleEnergy();
  G4double Eleak = Ebeam - calorEtot;
  runAct->fillPerEvent_3(calorEvis,calorEtot,Eleak);
  
  //nb of radiation lenght
  //
  runAct->fillNbRadLen(nbRadLen);  
  if (nbRadLen > 0.) histoManager->FillHisto(5,nbRadLen);
                 
  //parameters for trajectory visualisation
  //
  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                             ((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") trj->DrawTrajectory(100);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(100); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


