// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05CalorimeterSD.cc,v 1.2 1999-12-15 14:49:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN05CalorimeterSD.hh"
#include "ExN05CalorimeterHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

ExN05CalorimeterSD::ExN05CalorimeterSD(G4String name, G4int nCells, G4String colName)
:G4VSensitiveDetector(name),numberOfCells(nCells),HCID(-1)
{
  G4String HCname;
  collectionName.insert(HCname=colName);
  CellID = new int[numberOfCells];
}

ExN05CalorimeterSD::~ExN05CalorimeterSD()
{
  delete [] CellID;
}

void ExN05CalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  CalCollection = new ExN05CalorimeterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for(int j=0;j<numberOfCells;j++)
  {
    CellID[j] = -1;
  }
}

G4bool ExN05CalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep<=0.) return false;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int copyID = physVol->GetCopyNo();

  if(CellID[copyID]==-1)
  {
    ExN05CalorimeterHit* calHit = new ExN05CalorimeterHit(physVol->GetLogicalVolume());
    G4RotationMatrix rotM;
    if(physVol->GetObjectRotation()) rotM = *(physVol->GetObjectRotation());
    calHit->SetEdep( edep );
    calHit->SetPos( physVol->GetTranslation() );
    calHit->SetRot( rotM );
    int icell = CalCollection->insert( calHit );
    CellID[copyID] = icell - 1;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on CellID " << copyID << G4endl; }
  }
  else
  { 
    (*CalCollection)[CellID[copyID]]->AddEdep( edep );
    if(verboseLevel>0)
    { G4cout << " Energy added to CellID " << copyID << G4endl; }
  }

  return true;
}

void ExN05CalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, CalCollection );
}

void ExN05CalorimeterSD::clear()
{
} 

void ExN05CalorimeterSD::DrawAll()
{
} 

void ExN05CalorimeterSD::PrintAll()
{
} 




