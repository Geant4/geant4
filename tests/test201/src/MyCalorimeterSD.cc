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
// $Id: MyCalorimeterSD.cc,v 1.5 2006-06-29 21:47:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"

#include "MyCalorimeterSD.hh"
#include "MyCalorimeterHit.hh"

MyCalorimeterSD::MyCalorimeterSD(G4String name)
:G4VSensitiveDetector(name)
{
  collectionName.insert("CalCollection");
}

MyCalorimeterSD::~MyCalorimeterSD(){;}

void MyCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new MyCalorimeterHitsCollection(SensitiveDetectorName,
						  collectionName[0]); 
  for(int j=0;j<3;j++) {
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
    MyCalorimeterHit* calHit =
      new MyCalorimeterHit((G4VPhysicalVolume*)physVol);
    G4RotationMatrix rotM;
    if(physVol->GetObjectRotation()) rotM = *(physVol->GetObjectRotation());
    calHit->SetEdep( edep );
    calHit->SetPos( physVol->GetTranslation() );
    calHit->SetRot( rotM );
    int icell = CalCollection->insert(calHit);
    CellID[copyID] = icell;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on CellID " << copyID << G4endl; }
  }
  else
  { 
    //CalCollection->AddEdep( CellID[copyID], edep );
    if(verboseLevel>0)
    { G4cout << " Energy added to CellID " << copyID << G4endl; }
  }

  return true;
}

void MyCalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID, CalCollection );
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

