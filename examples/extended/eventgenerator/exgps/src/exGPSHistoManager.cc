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
/// \file eventgenerator/exgps/src/exGPSHistoManager.cc
/// \brief Implementation of the exGPSHistoManager class
//
//
// $Id: exGPSHistoManager.cc 74272 2013-10-02 14:48:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "exGPSHistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "exGPSHistoMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSHistoManager::exGPSHistoManager():
fMinpos(-10.),fMaxpos(10),fMineng(0.),fMaxeng(1000.),
fEnerHisto(0),fPosiXY(0),fPosiZX(0),fPosiYZ(0),fAnglCTP(0),fAnglTP(0)
{
  fFileName[0] = "expGPS";
  fFactoryOn = false;
  fMessenger = new exGPSHistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSHistoManager::~exGPSHistoManager()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSHistoManager::book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in exGPSHistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  G4String extension = analysisManager->GetFileType();
  fFileName[1] = fFileName[0] + "." + extension;
      
  // Create directories 
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");
    
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile(fFileName[0]);
  if (!fileOpen) {
    G4cout << "\n---> exGPSHistoManager::book(): cannot open " << fFileName[1]
           << G4endl;
    return;
  }
  
  // create selected histograms
  //
  analysisManager->SetFirstHistoId(1);
  analysisManager->SetFirstNtupleId(1);

  G4int id = analysisManager->CreateH1("1","Ekin primary",100,fMineng,fMaxeng);
  fEnerHisto= analysisManager->GetH1(id);

  id = analysisManager->CreateH2("2","Position XY",100,fMinpos,fMaxpos,
                                                   100,fMinpos,fMaxpos);
  fPosiXY= analysisManager->GetH2(id);

  id = analysisManager->CreateH2("3","Position YZ",100,fMinpos,fMaxpos,
                                                     100,fMinpos,fMaxpos);
  fPosiYZ= analysisManager->GetH2(id);

  id = analysisManager->CreateH2("4","Position ZX",100,fMinpos,fMaxpos,
                                                       100,fMinpos,fMaxpos);
  fPosiZX= analysisManager->GetH2(id);

  id =analysisManager->CreateH2("5","Source phi-std::cos(theta) distribution",
          360,0,360,100, -1, 1);
  fAnglCTP =analysisManager->GetH2(id);

  id =analysisManager->CreateH2("6","Source phi-theta distribution",
          360,0,360,180,0,180);
  fAnglTP = analysisManager->GetH2(id);

  // Create 1st ntuple (id = 1)
  //    
  analysisManager->CreateNtuple("MyTuple", "Primary Particle Tuple");
  id = analysisManager->CreateNtupleIColumn("particleID");
  id = analysisManager->CreateNtupleDColumn("Ekin");
  id = analysisManager->CreateNtupleDColumn("posX");
  id = analysisManager->CreateNtupleDColumn("posY");
  id = analysisManager->CreateNtupleDColumn("posZ");
  id = analysisManager->CreateNtupleDColumn("dirTheta");
  id = analysisManager->CreateNtupleDColumn("dirPhi");
  id = analysisManager->CreateNtupleDColumn("weight");
  analysisManager->FinishNtuple();
  
  fFactoryOn = true;       
  G4cout << "\n----> Histogram Tree is opened in " << fFileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSHistoManager::save()
{
  if (fFactoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();  
    G4cout << "\n----> Histogram Tree is saved in " << fFileName[1] << G4endl;
      
    delete G4AnalysisManager::Instance();
    fFactoryOn = false;
  }                    
}

void exGPSHistoManager::Fill(G4int PDGid, G4double e,
                                   G4double x, G4double y, G4double z,
                                   G4double t, G4double p, G4double w)
{
   fEnerHisto->fill(e, w);
   fEnerHisto->fill(e/MeV,w);
   fPosiXY->fill(x/cm,y/cm,w);
   fPosiZX->fill(z/cm,x/cm,w);
   fPosiYZ->fill(y/cm,z/cm,w);
   fAnglCTP->fill(p/deg,std::cos(t),w);
   fAnglTP->fill(p/deg,t/deg,w);

   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   analysisManager->FillNtupleIColumn(1,0,PDGid);
   analysisManager->FillNtupleDColumn(1,1,e);
   analysisManager->FillNtupleDColumn(1,2,x);
   analysisManager->FillNtupleDColumn(1,3,y);
   analysisManager->FillNtupleDColumn(1,4,z);
   analysisManager->FillNtupleDColumn(1,5,t);
   analysisManager->FillNtupleDColumn(1,6,p);
   analysisManager->FillNtupleDColumn(1,7,w);
   analysisManager->AddNtupleRow(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSHistoManager::PrintStatistic()
{
  if(fFactoryOn) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    /*G4cout
       << " EAbs : mean = " << G4BestUnit(fHistPt[1]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[1]->rms(),  "Energy") 
               << G4endl;
    G4cout                
       << " EGap : mean = " << G4BestUnit(fHistPt[2]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[2]->rms(),  "Energy") 
               << G4endl;
    G4cout 
       << " LAbs : mean = " << G4BestUnit(fHistPt[3]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[3]->rms(),  "Length") 
               << G4endl;
    G4cout 
       << " LGap : mean = " << G4BestUnit(fHistPt[4]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[4]->rms(),  "Length") 
               << G4endl;
     */
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


