// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterSD.cc,v 1.1 1999-04-16 10:32:36 johna Exp $
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
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="CalCollection");
}

MyCalorimeterSD::~MyCalorimeterSD(){;}

void MyCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  //CalCollection = new MyCalorimeterHitsCollection(collectionName[0],this); ??
  CalCollection = new MyCalorimeterHitsCollection(collectionName[0],""); 
  for(int j=0;j<3;j++)
  {
    CellID[j] = -1;
  }
}

G4bool MyCalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = 1.;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int copyID = physVol->GetCopyNo();

  if(CellID[copyID]==-1)
  {
    MyCalorimeterHit calHit((G4VPhysicalVolume*)physVol);
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
  // HCE->AddHitsCollection(CalCollection ); ??????????????????????
  HCE->AddHitsCollection(0, CalCollection );
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

