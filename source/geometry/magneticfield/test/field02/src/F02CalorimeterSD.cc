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
// $Id: F02CalorimeterSD.cc,v 1.2 2001-11-19 16:40:26 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F02CalorimeterSD.hh"

#include "F02CalorHit.hh"
#include "F02DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
  
#include "G4ios.hh"

/////////////////////////////////////////////////////////////////////////////////

F02CalorimeterSD::F02CalorimeterSD(G4String name,
                                   F02DetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("CalCollection");
  HitID = new G4int[500];
}

//////////////////////////////////////////////////////////////////////////////////

F02CalorimeterSD::~F02CalorimeterSD()
{
  delete [] HitID;
}

/////////////////////////////////////////////////////////////////////////////

void F02CalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  CalCollection = new F02CalorHitsCollection( SensitiveDetectorName,
                                              collectionName[0]      ); 
  for (G4int j=0;j<1; j++) 
  {
    HitID[j] = -1;
  }
}

////////////////////////////////////////////////////////////////////////////////

G4bool F02CalorimeterSD::ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  // if ( (aStep->GetTrack()->GetDefinition()->GetPDGCharge()) != 0.0 ) 
  {
    stepl = aStep->GetStepLength();
  }
  if ((edep == 0.) && (stepl == 0.) ) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();
 
  //theTouchable->MoveUpHistory();

  G4int F02Number = 0 ;

  if (HitID[F02Number]==-1)
  { 
    F02CalorHit* calHit = new F02CalorHit();

    if (physVol == Detector->GetAbsorber()) 
    {
      calHit->AddAbs(edep,stepl);
    }
    HitID[F02Number] = CalCollection->insert(calHit) - 1;

    if (verboseLevel>0)
    {  
      G4cout << " New Calorimeter Hit on F02: " << F02Number << G4endl;
    }
  }
  else
  { 
    if (physVol == Detector->GetAbsorber())
    {   
      (*CalCollection)[HitID[F02Number]]->AddAbs(edep,stepl);
    }
    if (verboseLevel>0)
    {    
      G4cout << " Energy added to F02: " << F02Number << G4endl;
    } 
  }    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02CalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02CalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void F02CalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

