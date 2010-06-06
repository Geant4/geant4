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
// $Id: EventAction.cc,v 1.9 2010-06-06 06:14:58 perl Exp $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det, RunAction* run,
                         PrimaryGeneratorAction* prim, HistoManager* hist)
:detector(det), runAct(run), primary(prim), histoManager(hist)
{ 
  trigger = false;
  Eseuil  = 10*keV;
  
  writeFile = false;
    
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
  G4int n1pxl   = detector->GetN1Pixels();
  G4int n2pxl   = detector->GetN2Pixels();
  G4int n1shift = detector->GetN1Shift();
        
  // code for trigger conditions :
  // 1 and only 1 pixel fired per layer
  //
  if (trigger) {
    for (G4int i1=0; i1<n1pxl; i1++) {
      //count number of pixels fired
      G4int count = 0;  
      for (G4int i2=0; i2<n2pxl; i2++) {
        G4int k = i1*n1shift + i2;
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
  for (G4int i1=0; i1<n1pxl; i1++) {
    //sum energy per readout layer  
    G4double layerEvis = 0.;
    G4double layerEtot = 0.;  
    for (G4int i2=0; i2<n2pxl; i2++) {
      G4int k = i1*n1shift + i2;
      runAct->fillPerEvent_1(k,visibleEnergy[k],totalEnergy[k]);      
      layerEvis += visibleEnergy[k];
      layerEtot += totalEnergy[k];
      calorEvis += visibleEnergy[k];
      calorEtot += totalEnergy[k];		      
    }      
    runAct->fillPerEvent_2(i1,layerEvis,layerEtot);
    if (layerEvis > 0.) histoManager->FillNtuple(1, i1, layerEvis);
    if (layerEtot > 0.) histoManager->FillNtuple(1, n1pxl+i1, layerEtot);
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
  
  //write file of pixels
  //
  if (writeFile) WritePixels(evt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SetWriteFile(G4bool val)    
{
  writeFile = val;
  runAct->SetWriteFile(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
#include <fstream>
#include "G4RunManager.hh"
#include "G4Run.hh"

void EventAction::WritePixels(const G4Event* evt)
{
  // event is appended onto file created at BeginOfRun
  //
  G4String name = histoManager->GetFileName(); 
  G4String fileName = name + ".pixels.ascii";

  std::ofstream File(fileName, std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
    
  //write event number  
  //
  File << evt->GetEventID() << G4endl;
  
  //gun particle informations
  //
  G4ParticleGun* gun = primary->GetParticleGun();
  G4double ekin = gun->GetParticleEnergy();
  G4ThreeVector direction = gun->GetParticleMomentumDirection();
  G4ThreeVector position  = gun->GetParticlePosition();
  File << ekin << " " << direction << " " << position << G4endl;  
    
  //count nb of fired pixels
  //
  G4int firedPixels = 0;
  G4int nbOfPixels = detector->GetSizeVectorPixels();  
  for (G4int k=0; k<nbOfPixels; k++) {
    if (totalEnergy[k] > 0.0) firedPixels++;
  }         
  File << firedPixels << G4endl;
  
  //write pixels
  //
  for (G4int k=0; k<nbOfPixels; k++) {
    if (totalEnergy[k] > 0.0) 
    File << k << " " << visibleEnergy[k] << " " << totalEnergy[k] << " "; 
  }            
  File << G4endl;
    
  // restaure default formats
  File.setf(mode,std::ios::floatfield);
  File.precision(prec);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


