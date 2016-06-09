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
// $Id: RunAction.cc,v 1.1 2003/08/11 10:21:32 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(HistoManager* histo)
:histoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //initialisation
  EnergyDeposit  = EnergyDeposit2  = 0.;
  TrakLenCharged = TrakLenCharged2 = 0.;
  TrakLenNeutral = TrakLenNeutral2 = 0.;
  nbStepsCharged = nbStepsCharged2 = 0.;
  nbStepsNeutral = nbStepsNeutral2 = 0.;
  
  nbGamma = nbElect = nbPosit = 0;
  
  Transmit[0] = Transmit[1] = Reflect[0] = Reflect[1] = 0;
          
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  // compute mean and rms
  //
  G4int TotNbofEvents = aRun->GetNumberOfEvent();
  if (TotNbofEvents == 0) return;
  
  EnergyDeposit /= TotNbofEvents;
  G4double rmsEdep = EnergyDeposit2/TotNbofEvents - EnergyDeposit*EnergyDeposit;
  if (rmsEdep>0.) rmsEdep = sqrt(rmsEdep/TotNbofEvents);
  else            rmsEdep = 0.;
  
  TrakLenCharged /= TotNbofEvents;
  G4double rmsTLCh = TrakLenCharged2/TotNbofEvents 
                                                - TrakLenCharged*TrakLenCharged;
  if (rmsTLCh>0.) rmsTLCh = sqrt(rmsTLCh/TotNbofEvents);
  else            rmsTLCh = 0.;
 
  TrakLenNeutral /= TotNbofEvents;
  G4double rmsTLNe = TrakLenNeutral2/TotNbofEvents 
                                                - TrakLenNeutral*TrakLenNeutral;
  if (rmsTLNe>0.) rmsTLNe = sqrt(rmsTLNe/TotNbofEvents);
  else            rmsTLNe = 0.;
  
  nbStepsCharged /= TotNbofEvents;
  G4double rmsStCh = nbStepsCharged2/TotNbofEvents 
                                                - nbStepsCharged*nbStepsCharged;
  if (rmsStCh>0.) rmsStCh = sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStCh = 0.;  
  
  nbStepsNeutral /= TotNbofEvents;
  G4double rmsStNe = nbStepsNeutral2/TotNbofEvents 
                                                - nbStepsNeutral*nbStepsNeutral;
  if (rmsStNe>0.) rmsStNe = sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStNe = 0.;
  
  G4double Gamma = (double)nbGamma/TotNbofEvents; 
  G4double Elect = (double)nbElect/TotNbofEvents;
  G4double Posit = (double)nbPosit/TotNbofEvents;
  
  G4double transmit[2];
  transmit[0] = 100.*Transmit[0]/TotNbofEvents;
  transmit[1] = 100.*Transmit[1]/TotNbofEvents;
    
  G4double reflect[2];
  reflect[0] = 100.*Reflect[0]/TotNbofEvents;
  reflect[1] = 100.*Reflect[1]/TotNbofEvents;
       
 G4cout << "\n ======================== run summary ======================\n"; 

 G4int prec = G4cout.precision(4);
 
 G4cout << "\n Number of Events = " << TotNbofEvents << G4endl;
  
 G4cout << "\n Total energy deposit in absorber per event = " 
        << G4BestUnit(EnergyDeposit,"Energy") << " +- "
        << G4BestUnit(rmsEdep,      "Energy") << G4endl;
	
 G4cout << "\n Total track length (charged) in absorber per event = " 
        << G4BestUnit(TrakLenCharged,"Length") << " +- "
        << G4BestUnit(rmsTLCh,       "Length") << G4endl;
	
 G4cout << " Total track length (neutral) in absorber per event = " 
        << G4BestUnit(TrakLenNeutral,"Length") << " +- "
        << G4BestUnit(rmsTLNe,       "Length") << G4endl;
			   
 G4cout << "\n Number of steps (charged) in absorber per event = " 
        << nbStepsCharged << " +- " << rmsStCh << G4endl;
	
 G4cout << " Number of steps (neutral) in absorber per event = " 
        << nbStepsNeutral << " +- " << rmsStNe << G4endl;
	
 G4cout << "\n Number of secondaries per event : Gammas = " << Gamma 
        << ";   electrons = " << Elect 
	<< ";   positrons = " << Posit << G4endl;
	
 G4cout << "\n Number of events with the primary particle transmitted = " 
        << transmit[1] << " %" << G4endl;
	
 G4cout << " Number of events with at least  1 particle transmitted " 
        << "(same charge as primary) = " << transmit[0] << " %" << G4endl;
	
 G4cout << "\n Number of events with the primary particle reflected = " 
        << reflect[1] << " %" << G4endl;
	
 G4cout << " Number of events with at least  1 particle reflected " 
        << "(same charge as primary) = " << reflect[0] << " %" << G4endl;
					
 G4cout.precision(prec);
  
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
   
  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
