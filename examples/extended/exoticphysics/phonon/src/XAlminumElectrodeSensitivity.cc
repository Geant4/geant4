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
/// \file exoticphysics/phonon/src/XAlminumElectrodeSensitivity.cc
/// \brief Implementation of the XAlminumElectrodeSensitivity class
//
// $Id$
//
#include "XAlminumElectrodeSensitivity.hh"

#include "XAlminumElectrodeHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


using namespace std;

XAlminumElectrodeHitsCollection* XAlminumElectrodeSensitivity::hitsCollection = NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAlminumElectrodeSensitivity::XAlminumElectrodeSensitivity(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="XAlminumElectrodeHit");
  HCID = -1;
  fWriter.open("caustic.ssv", fstream::out | fstream::ate);
  fWriter2.open("timing.ssv", fstream::out | fstream::ate);

  if(!fWriter.is_open()){
    G4cout<<"\nXAlminumElectrodeSensitivity::Constructor:";
    G4cout<<"\n\tFailed to open caustic.ssv for appending data.";
    G4cout<<"\n\tCreating caustic.ssv";
    fWriter.open("caustic.ssv");
  }

  if(!fWriter2.is_open()){
    G4cout<<"\nXAlminumElectrodeSensitivity::Constructor: ";
    G4cout<<"\n\tFailed to open timing.ssv for appending data.";
    G4cout<<"\n\tCreating timing.ssv.";
    fWriter2.open("timing.ssv");
  }

  if(!(fWriter.is_open() && fWriter2.is_open())){
    G4cout<<"\nXAlminumElectrodeSensitivity::Constructor: ERROR: COULD NOT CREATE OUTPUT FILES FOR WRITING";
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAlminumElectrodeSensitivity::~XAlminumElectrodeSensitivity(){
  fWriter.close();
  fWriter2.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAlminumElectrodeHitsCollection* XAlminumElectrodeSensitivity::GetHitsCollection(){
  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XAlminumElectrodeSensitivity::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new XAlminumElectrodeHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XAlminumElectrodeSensitivity::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
  //if(aStep->GetTrack()->GetDefinition()!=Phonon::PhononDefinition()) return true;
  G4double edp = aStep->GetNonIonizingEnergyDeposit();
  if(edp==0.) return true;

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4ThreeVector fWorldPos = postStepPoint->GetPosition();
  G4ThreeVector fLocalPos
    = theTouchable->GetHistory()->GetTopTransform().TransformPoint(fWorldPos);

  XAlminumElectrodeHit* aHit = new XAlminumElectrodeHit();
  aHit->SetTime(postStepPoint->GetGlobalTime());
  aHit->SetEDep(edp);
  aHit->SetWorldPos(fWorldPos);
  aHit->SetLocalPos(fLocalPos);

  hitsCollection->insert(aHit);

  fWriter<<"\n"<<fWorldPos.getX()/mm
         <<","<<fWorldPos.getY()/mm
         <<","<<fWorldPos.getZ()/mm;

  fWriter2<<"\n"<<postStepPoint->GetGlobalTime()/ns<<" "
          <<aHit->GetEDep()/eV;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XAlminumElectrodeSensitivity::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{;}

