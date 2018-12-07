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
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
// 15 Jul 2003 Alfonso Mantero        "DetectorType" use integration
// 23 Sep 2003 Alfonso Mantero        differnt geometries integration
//
// -------------------------------------------------------------------

#include "XrayFluoEventAction.hh"
#include "XrayFluoSensorHit.hh"
#include "XrayFluoEventActionMessenger.hh"
#include "XrayFluoDataSet.hh"
#include "XrayFluoAnalysisManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoEventAction::XrayFluoEventAction(const XrayFluoDetectorConstruction* det)
  :drawFlag("all"),
   HPGeCollID(0),
   eventMessenger(0),
   printModulo(1),
   detectorType(0)
{
  eventMessenger = new XrayFluoEventActionMessenger(this);
  
  if (!(det->GetPhaseSpaceFlag()) ){
    detectorType = det->GetDetectorType();
    HPGeCollID=-1;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoEventAction::XrayFluoEventAction(const XrayFluoPlaneDetectorConstruction* det)
  :drawFlag("all"),
   HPGeCollID(-1),
   eventMessenger(0),
   printModulo(1),
   detectorType(0)
{
  eventMessenger = new XrayFluoEventActionMessenger(this);
  detectorType = det->GetDetectorType();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoEventAction::XrayFluoEventAction(const XrayFluoMercuryDetectorConstruction* det)
  :drawFlag("all"),
   HPGeCollID(-1),
   eventMessenger(0),
   printModulo(1),
   detectorType(0)
{
  eventMessenger = new XrayFluoEventActionMessenger(this);
  detectorType = det->GetDetectorType();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoEventAction::~XrayFluoEventAction()
{
   delete eventMessenger;
   eventMessenger = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoEventAction::BeginOfEventAction(const G4Event* evt)
{

  G4int eventNumber = (evt->GetEventID())+1;
  if ( eventNumber == 1){

  G4cout << "# = 100000 events" << G4endl;
  G4cout << "1--------+---------+---------+---------+---------5e6"<<G4endl;
  }

  if ( ((eventNumber) % 100000) == 0 )  {

    if ( eventNumber % (G4int)5e6 != 0 ) G4cout << "#" << std::flush;
    else G4cout << "#"<< G4endl;   
  }

  if (HPGeCollID==-1)    
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      HPGeCollID = SDman->GetCollectionID("HPGeCollection");
      //the pointer points to the ID number of the sensitive detector
    }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoEventAction::EndOfEventAction(const G4Event* evt)
{
  
  if (detectorType) {
    
    // extracted from hits, compute the total energy deposit (and total charged
    // track length) 
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    
    XrayFluoSensorHitsCollection* HPGeHC = 0;
    G4int n_hit = 0;
    G4double totEnergyDetect=0., totEnergy=0., energyD=0.;
    
    if (HCE) HPGeHC = 
      (XrayFluoSensorHitsCollection*)(HCE->GetHC(HPGeCollID));

    if(HPGeHC)            
      {
	n_hit = HPGeHC->entries();
    
	for (G4int i=0;i<n_hit;i++)
	  {
	    totEnergy += (*HPGeHC)[i]->GetEdepTot(); 
	    energyD = detectorType->ResponseFunction(totEnergy);	    
	    XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
	    analysis->analyseEnergyDep(energyD);
	    totEnergyDetect += energyD;	    	    
	  }
      }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XrayFluoEventAction::RandomCut(G4double energy)
  
{
  G4double efficiency = 1.;
  G4double F = 0.15;
  G4double epsilon = 2.96 * eV;
  G4double deltaE = 220 * eV;
  G4double EdepDetect = 0.;

  G4double  Random = G4UniformRand(); 

  if ( Random<efficiency )
    {
      G4double sigma = std::sqrt(F*epsilon*energy+std::pow(deltaE/2355,2));      
      EdepDetect = G4RandGauss::shoot(energy, sigma );
      
    }
  else {EdepDetect = 0.;}

  return   EdepDetect;  
}


