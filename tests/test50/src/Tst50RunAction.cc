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
// $Id: Tst50RunAction.cc,v 1.6 2003-01-17 17:14:15 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst50AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
Tst50RunAction::Tst50RunAction()
{
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50RunAction::~Tst50RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl; 

  runID=aRun->GetRunID();


#ifdef G4ANALYSIS_USE   
Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
   analysis->book();
#endif

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
  number=0;
  numberTransp=0;
  numberRay=0;
  numberPh=0; 
  numberCo=0;
  numberPair=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int Tst50RunAction::GetRun_ID()
{
  return runID;
}
void Tst50RunAction::EndOfRunAction(const G4Run*)
{
#ifdef G4ANALYSIS_USE
 Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();

  analysis->finish();
#endif
  if (G4VVisManager::GetConcreteInstance())
    {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }
  
  //G4cout<<"---------- Particelle primarie trasmesse -------------- "<<number<<G4endl;
  if (numberTransp==0 && numberRay==0 && numberPh==0 && numberCo==0 &&numberPair==0){;}
else
    { 
G4cout<<"---------- Particelle primarie gamma -------------- "<<number<<G4endl;
 G4cout<<"--------- Processi verificatesi----------------"<<G4endl;
  G4cout<<"Number of transmitted gamma: "<<number<<G4endl;
  G4cout<<numberTransp <<" processo di trasporto"<< G4endl;
  G4cout<<numberRay<<" processi Rayleigh"<<G4endl;
  G4cout<<numberPh<< " processi fotoelettrici"<<G4endl;
  G4cout<<numberCo<< " processi Compton"<< G4endl;
  G4cout<<numberPair<< " processi di produzione di coppie"<< G4endl;}
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  Tst50RunAction::Trans_number()
{
  number= number+1;
}

void  Tst50RunAction::primary_processes(G4int i)
{
  if( i==1) numberTransp=numberTransp+1;
  if( i==2) numberRay= numberRay+1;
  if( i==3) numberPh= numberPh+1;
  if( i==4) numberCo= numberCo+1;
  if( i==5) numberPair= numberPair+1;
}
