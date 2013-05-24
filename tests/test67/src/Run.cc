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
// $Id: Run.cc 66241 2012-12-13 18:34:42Z gunter $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Run.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "Hits.hh"
#include "G4SDManager.hh"

#ifdef G4_USE_ROOT
#include "ROOTAnalysis.hh"
#endif

#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run() :
  counter(0),counterTot(0),fPrimaryEnergy(-1),
  hcID(-1)
{
  //Retrieve primary energy, if available
  const PrimaryGeneratorAction* primGen = 
    dynamic_cast<const PrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  //Primary generator not yet built -> the instance was created by 
  //the master RunAction in MT mode.
  if (!primGen)
    return;
 
  fPrimaryEnergy = primGen->GetPrimaryEnergy();
  
  //G4cout << "Primary energy: " << fPrimaryEnergy/keV << " keV " << G4endl;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* evt)
{
  //First, get the total event energy
  G4HCofThisEvent *HCE=evt->GetHCofThisEvent(); 

  if (hcID <0)
    {
      G4SDManager *SDman = G4SDManager::GetSDMpointer();
      hcID = SDman->GetCollectionID("HitsCollection");      
    }

  HitsCollection *CHC = 0; 
  size_t n_hit = 0;
  G4double etot = 0.;
 
  if (HCE) 
    CHC=(HitsCollection*) (HCE->GetHC(hcID));

  if (CHC) 
    {
      n_hit=CHC->entries();
      for (size_t i=0;i<n_hit;i++)	
	etot += (*CHC)[i]->GetEdep(); //le entries sono indicizzate
         
    }

  //Second, assess if the counters have to tick
  if (etot > 0.1*eV) 
    {
      counterTot++;

#ifdef G4_USE_ROOT
      ROOTAnalysis::getInstance()->AddEventEnergy(GetRunID(),etot);
#endif

      if (std::fabs(GetPrimaryEnergy() - etot) < 0.001*GetPrimaryEnergy()) 
	//within 0.1%
	counter++;
    }

  //Third, close it
  G4Run::RecordEvent(evt);

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void Run::Merge(const G4Run* aRun)
{
  const Run* localRun = static_cast<const Run*>(aRun);

  //increment counters
  counter += localRun->counter;
  counterTot += localRun->counterTot;
  fPrimaryEnergy = localRun->fPrimaryEnergy; //copy from any of the slaves

  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

G4double Run::GetPrimaryEnergy() const
{
  return fPrimaryEnergy;
}

