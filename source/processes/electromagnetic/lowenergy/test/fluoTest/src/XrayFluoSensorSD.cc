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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoSensorSD.hh"
#include "XrayFluoSensorHit.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorSD::XrayFluoSensorSD(G4String name,
                                   XrayFluoDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("SenCollection");
  HitID = new G4int[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSensorSD::~XrayFluoSensorSD()
{
  delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorSD::Initialize(G4HCofThisEvent*HCE)
  
{
  SenCollection = new XrayFluoSensorHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  
  {HitID[0] = -1;};
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XrayFluoSensorSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
  
  if ((edep==0.)&&(stepl==0.)) return false; 
  
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  //theTouchable->MoveUpHistory();
  
  if (HitID[0]==-1)
    { 
      XrayFluoSensorHit* senHit = new XrayFluoSensorHit();
      if (physVol == Detector->GetSensor()) 
	senHit->AddSen(edep,stepl);
      HitID[0] = SenCollection->insert(senHit) - 1;
      if (verboseLevel>0)
        G4cout << " New Sensor Hit  " << G4endl;
    }
  
  else
    { 
      if (physVol == Detector->GetSensor())
	{(*SenCollection)[HitID[0]]->AddSen(edep,stepl);}
      if (verboseLevel>0)
        {G4cout << " Energy added to the detector "  << G4endl;} 
    }
  
  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
  HCE->AddHitsCollection(HCID,SenCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSensorSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



