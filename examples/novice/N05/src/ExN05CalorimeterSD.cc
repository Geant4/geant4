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
// $Id: ExN05CalorimeterSD.cc,v 1.6 2006/06/29 17:53:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

#include "ExN05CalorimeterSD.hh"
#include "ExN05CalorimeterHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

ExN05CalorimeterSD::ExN05CalorimeterSD( G4String name,
                                        G4int nCells,
                                        G4String colName )
  : G4VSensitiveDetector(name),
    numberOfCells(nCells),
    HCID(-1)
{
  G4String HCname;
  collectionName.insert(HCname=colName);
  CellID = new G4int[numberOfCells];
}

ExN05CalorimeterSD::~ExN05CalorimeterSD()
{
  delete [] CellID;
}

void ExN05CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new ExN05CalorimeterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for(G4int j=0;j<numberOfCells;j++)
  {
    CellID[j] = -1;
  }
}

G4bool ExN05CalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep<=0.) return false;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4int copyID = physVol->GetCopyNo();

  if(CellID[copyID]==-1)
  {
    ExN05CalorimeterHit* calHit =
      new ExN05CalorimeterHit(physVol->GetLogicalVolume());
    G4RotationMatrix rotM;
    if(physVol->GetObjectRotation()) rotM = *(physVol->GetObjectRotation());
    calHit->SetEdep( edep );
    calHit->SetPos( physVol->GetTranslation() );
    calHit->SetRot( rotM );
    G4int icell = CalCollection->insert( calHit );
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
