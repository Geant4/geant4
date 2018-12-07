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
/// \file Par01/src/Par01CalorimeterSD.cc
/// \brief Implementation of the Par01CalorimeterSD class
//
//
//

#include "Par01CalorimeterSD.hh"
#include "Par01CalorimeterHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01CalorimeterSD::Par01CalorimeterSD( G4String name,
                                        G4int nCells,
                                        G4String colName )
  : G4VSensitiveDetector(name),
    fNumberOfCells(nCells),
    fHCID(-1)
{
  G4String HCname;
  collectionName.insert(HCname=colName);
  fCellID = new G4int[fNumberOfCells];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01CalorimeterSD::~Par01CalorimeterSD()
{
  delete [] fCellID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  fCalCollection = new Par01CalorimeterHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  for(G4int j=0;j<fNumberOfCells;j++)
    {
      fCellID[j] = -1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par01CalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep<=0.) return false;
  
  G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  const G4VPhysicalVolume* physVol = hist->GetVolume();
  G4int copyID = hist->GetReplicaNumber();
  
  if(fCellID[copyID]==-1)
    {
      Par01CalorimeterHit* calHit =
        new Par01CalorimeterHit(physVol->GetLogicalVolume());
      calHit->SetEdep( edep );
      G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
      aTrans.Invert();
      calHit->SetPos(aTrans.NetTranslation());
      calHit->SetRot(aTrans.NetRotation());
      G4int icell = fCalCollection->insert( calHit );
      fCellID[copyID] = icell - 1;
      if(verboseLevel>0)
        { G4cout << " New Calorimeter Hit on CellID " << copyID << G4endl; }
    }
  else
    { 
      (*fCalCollection)[fCellID[copyID]]->AddEdep( edep );
      if(verboseLevel>0)
        { G4cout << " Energy added to CellID " << copyID << G4endl; }
    }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01CalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  if(fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( fHCID, fCalCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01CalorimeterSD::clear()
{
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01CalorimeterSD::DrawAll()
{
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01CalorimeterSD::PrintAll()
{
} 
