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
// $Id$
//
/// \file analysis/A01/src/A01HadCalorimeter.cc
/// \brief Implementation of the A01HadCalorimeter class

#include "A01HadCalorimeter.hh"
#include "A01HadCalorimeterHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

A01HadCalorimeter::A01HadCalorimeter(G4String name)
:G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    G4String HCname = "HadCalorimeterColl";
    collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

A01HadCalorimeter::~A01HadCalorimeter(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void A01HadCalorimeter::Initialize(G4HCofThisEvent*HCE)
{
    fHitsCollection = new A01HadCalorimeterHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if(fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    HCE->AddHitsCollection(fHCID,fHitsCollection);
    
    // fill calorimeter hits with zero energy deposition
    for(int iColumn=0;iColumn<10;iColumn++)
        for(int iRow=0;iRow<2;iRow++)
        {
            A01HadCalorimeterHit* aHit = new A01HadCalorimeterHit();
            fHitsCollection->insert( aHit );
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool A01HadCalorimeter::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    if(edep==0.) return true;
    
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* theCellPhysical = theTouchable->GetVolume(2);
    G4int rowNo = theCellPhysical->GetCopyNo();
    G4VPhysicalVolume* theColumnPhysical = theTouchable->GetVolume(3);
    G4int columnNo = theColumnPhysical->GetCopyNo();
    G4int hitID = 2*columnNo+rowNo;
    A01HadCalorimeterHit* aHit = (*fHitsCollection)[hitID];
    
    // check if it is first touch
    if(aHit->GetColumnID()<0)
    {
        aHit->SetColumnID(columnNo);
        aHit->SetRowID(rowNo);
        G4int depth = theTouchable->GetHistory()->GetDepth();
        G4AffineTransform aTrans = theTouchable->GetHistory()->GetTransform(depth-2);
        aTrans.Invert();
        aHit->SetRot(aTrans.NetRotation());
        aHit->SetPos(aTrans.NetTranslation());
    }
    // add energy deposition
    aHit->AddEdep(edep);
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
