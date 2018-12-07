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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelCalorimeterSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
#include "G4RunManager.hh"
#include "GammaRayTelCalorimeterSD.hh"
#include "GammaRayTelCalorimeterHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterSD::GammaRayTelCalorimeterSD(G4String name):G4VSensitiveDetector(name)
{
 G4RunManager* runManager = G4RunManager::GetRunManager();
  Detector =
    (GammaRayTelDetectorConstruction*)(runManager->GetUserDetectorConstruction());
  
  NbOfCALBars  = Detector->GetNbOfCALBars();
  NbOfCALLayers  = Detector->GetNbOfCALLayers();

  //G4cout <<  NbOfCALBars << " bars " << G4endl;
  //G4cout <<  NbOfCALLayers << " layers " << G4endl;
  
  NbOfCALChannels = NbOfCALBars*NbOfCALLayers;
  
  ChitXID = new G4int[NbOfCALChannels];
  ChitYID = new G4int[NbOfCALChannels];
  collectionName.insert("CalorimeterCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterSD::~GammaRayTelCalorimeterSD()
{
  delete [] ChitXID;
  delete [] ChitYID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  CalorimeterCollection = new GammaRayTelCalorimeterHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  for (G4int i=0;i<NbOfCALChannels;i++)
      {
	ChitXID[i] = -1;
	ChitYID[i] = -1;
      };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelCalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{ 
  
  G4double edep = 0.;
  edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;      
  
  // This TouchableHistory is used to obtain the physical volume
  // of the hit
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  
  G4VPhysicalVolume* cal_bar = theTouchable->GetVolume();  
  G4VPhysicalVolume* cal_plane = theTouchable->GetVolume(1);

  G4int CALBarNumber=cal_bar->GetCopyNo();
  G4String CALBarName = cal_bar->GetName();
  
  G4int PlaneNumber = 0;
  PlaneNumber=cal_plane->GetCopyNo();
  G4String PlaneName = cal_plane->GetName();


  G4int NChannel = 0;
  
  NChannel = PlaneNumber * NbOfCALBars + CALBarNumber; 
  
  if (PlaneName == "CALLayerX" )
    
    // The hit is on an X CsI plane
    
    {
      // This is a new hit
      if (ChitXID[NChannel]==-1)
	{       
	  GammaRayTelCalorimeterHit* CalorimeterHit = new GammaRayTelCalorimeterHit;
	  CalorimeterHit->SetCALType(1);
	  CalorimeterHit->AddEnergy(edep);
	  CalorimeterHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  CalorimeterHit->SetCALPlaneNumber(PlaneNumber);
	  CalorimeterHit->SetCALBarNumber(CALBarNumber);
	  ChitXID[NChannel] = 
	    CalorimeterCollection->insert(CalorimeterHit) -1;
	}
      else // This is not new
	{
	  (*CalorimeterCollection)
	    [ChitXID[NChannel]]->AddEnergy(edep);
	}
    }
 
  if (PlaneName == "CALLayerY")
    // The hit is on an Y CsI plane    
    {   
      // This is a new hit
      if (ChitYID[NChannel]==-1)
	{       
	  GammaRayTelCalorimeterHit* CalorimeterHit 
	    = new GammaRayTelCalorimeterHit;
	  CalorimeterHit->SetCALType(0);
	  CalorimeterHit->AddEnergy(edep);
	  CalorimeterHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  CalorimeterHit->SetCALPlaneNumber(PlaneNumber);
	  CalorimeterHit->SetCALBarNumber(CALBarNumber);
	  ChitYID[NChannel] = 
	    CalorimeterCollection->insert(CalorimeterHit)-1;
	}
      else // This is not new
	{
	  (*CalorimeterCollection)
	    [ChitYID[NChannel]]->AddEnergy(edep);
	}
    }
  
  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void GammaRayTelCalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,CalorimeterCollection);


  for (G4int i=0;i<NbOfCALChannels;i++) 
    {
      ChitXID[i] = -1;
      ChitYID[i] = -1;
    };
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














