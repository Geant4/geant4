// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterSD.cc,v 1.1 1999-01-07 16:04:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyCalorimeterSD.hh"
#include "MyCalorimeterHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"

MyCalorimeterSD::MyCalorimeterSD(G4String name)
:G4VSensitiveDetector(name),numberOfCells(2550)
{
  G4String HCname;
  collectionName.insert(HCname="CalCollection");
  CellID = new int[numberOfCells];
}

MyCalorimeterSD::~MyCalorimeterSD()
{
  delete [] CellID;
}

void MyCalorimeterSD::Initialize()
{
  CalCollection = new MyCalorimeterHitsCollection(collectionName[0],this); 
  for(int j=0;j<numberOfCells;j++)
  {
    CellID[j] = -1;
  }
}

G4bool MyCalorimeterSD::ProcessHits(G4Step*aStep)
{
  G4double edep = aStep->GetDeltaEnergy();

  if(edep<=0.) return false;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int copyID = physVol->GetCopyNo();

  if(CellID[copyID]==-1)
  {
    MyCalorimeterHit calHit(physVol->GetLogicalVolume());
    G4RotationMatrix rotM;
    if(physVol->GetObjectRotation()) rotM = *(physVol->GetObjectRotation());
    calHit.SetEdep( edep );
    calHit.SetPos( physVol->GetTranslation() );
    calHit.SetRot( rotM );
    int icell = CalCollection->insert( &calHit );
    CellID[copyID] = icell;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on CellID " << copyID << endl; }
  }
  else
  { 
    CalCollection->AddEdep( CellID[copyID], edep );
    if(verboseLevel>0)
    { G4cout << " Energy added to CellID " << copyID << endl; }
  }

  return true;
}

void MyCalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  HCE->AddHitsCollection( CalCollection );
}

void MyCalorimeterSD::clear()
{
} 

void MyCalorimeterSD::DrawAll()
{
} 

void MyCalorimeterSD::PrintAll()
{
} 

