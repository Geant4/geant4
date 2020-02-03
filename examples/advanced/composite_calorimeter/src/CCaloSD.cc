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
///////////////////////////////////////////////////////////////////////////////
// File: CCaloSD.cc
// Description: Stores hits of calorimetric type in appropriate container
///////////////////////////////////////////////////////////////////////////////

#include "CCaloSD.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "G4SDManager.hh"
#include "G4VTouchable.hh"
#include "CCalVOrganization.hh"
#include "CCalSDList.hh"

//#define debug
//#define ddebug
 
CCaloSD::CCaloSD(G4String name, CCalVOrganization* numberingScheme):
  G4VSensitiveDetector(name), HCID(-1), SDname(name), theHC(0),
  CurrentHit(0), theTrack(0), CurrentPV(0), PreviousPV(0), UnitID(0), 
  PreviousUnitID(0), PreStepPoint(0), PostStepPoint(0), 
  theDescription(numberingScheme) {
  
  collectionName.insert(name);
  
  G4cout << "*******************************************************" << G4endl;
  G4cout << "*                                                     *" << G4endl;
  G4cout << "* Constructing a CCaloSD  with name " << name            << G4endl;
  G4cout << "*                                                     *" << G4endl;
  G4cout << "*******************************************************" << G4endl;

  CCalSDList::getInstance()->addCalo(name);
}

CCaloSD::~CCaloSD() {
  delete theDescription;
}

void CCaloSD::Initialize(G4HCofThisEvent*HCE) {

#ifdef debug
  G4cout << "CCaloSD : Initialize called for " << SDname << G4endl;
#endif
  //This initialization is performed at the beginning of an event
  //------------------------------------------------------------

  theHC = new CCalG4HitCollection(SDname, collectionName[0]);
  if (HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, theHC );

  TSID = -2;
  PrimID = -2;
  /////PrimaryID = -99;  <--- initialized by StackingAction.
}

G4bool CCaloSD::ProcessHits(G4Step*aStep, G4TouchableHistory*) {

  if(aStep->GetTotalEnergyDeposit() > CLHEP::eV) {
    getStepInfo(aStep);
    if (hitExists() == false && EdepositEM+EdepositEHAD>0.f) {
      createNewHit();
    }
  }
  return true;
} 

void CCaloSD::getStepInfo(const G4Step* aStep) {
  
  PreStepPoint = aStep->GetPreStepPoint(); 
  PostStepPoint= aStep->GetPostStepPoint(); 
  HitPoint     = PreStepPoint->GetPosition();        

  theTrack = aStep->GetTrack();   
  CurrentPV= PreStepPoint->GetPhysicalVolume();

  G4String pname = theTrack->GetDefinition()->GetParticleName();
  G4double de = aStep->GetTotalEnergyDeposit();
  //G4cout << "##### Process new step dE= " << de
  //         << "  " << pname << " inside " << GetName() << G4endl;

  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  G4double weight = 1.;
  if (touch->GetVolume(0)->GetName() == "CrystalMatrixCrystal") {
    weight = curve_LY(PreStepPoint);
  }

  if (pname == "e-" || pname == "e+" || pname == "gamma" ){
    EdepositEM   = weight*de;
    EdepositEHAD = 0.f;
  } else {
    EdepositEM   = 0.f;
    EdepositEHAD = weight*de;
  }

  TSlice = (PostStepPoint) ? PostStepPoint->GetGlobalTime()/nanosecond : 0.;
  //G4cout << "     W= " << weight << " T= " << TSlice << G4endl;
  G4int it = (G4int)TSlice;
  TSliceID = std::min(it, 10000);
  //G4cout << " tID= " <<  TSliceID << G4endl;

  UnitID = (theDescription) ? theDescription->GetUnitID(aStep) : 0;
}

G4bool CCaloSD::hitExists() {
   
  if ( CurrentPV==PreviousPV && PrimaryID == PrimID && TSliceID == TSID &&
       UnitID==PreviousUnitID) {
    updateHit();
    return true;
  }
   
  if (PrimaryID != PrimID) { ResetForNewPrimary(); }
   
  // look in HC whether a hit with the same primID,UnitID,TSliceID 
  // exists already:
   
  G4bool found = false;
  for (std::size_t j=0; j<theHC->entries(); ++j) {

    CCalG4Hit* aPreviousHit = (*theHC)[j];
    if (aPreviousHit->getTrackID()  == PrimaryID &&
        aPreviousHit->getTimeSliceID() == TSliceID  &&
        aPreviousHit->getUnitID()== UnitID       ) {
      CurrentHit = aPreviousHit;
      found = true;
      break;
    }
  }          

  if (found) {
    updateHit();
    return true;
  } else {
    return false;
  }    
}

void CCaloSD::ResetForNewPrimary() {
  
  EntrancePoint = SetToLocal(HitPoint);
  IncidentEnergy = PreStepPoint->GetKineticEnergy();

}

void CCaloSD::StoreHit(CCalG4Hit* hit){

  if (PrimID<0) return;
  if (hit == 0 ) {
    G4cout << "CCaloSD: hit to be stored is NULL !!" <<G4endl;
    return;
  }

  theHC->insert( hit );
}

void CCaloSD::createNewHit() {

#ifdef debug
  G4int currentCopyNo = -999;
  G4int motherCopyNo  = -999;
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)( theTrack->GetTouchable() );
  if ( theTouchable ) {
    currentCopyNo = theTouchable->GetReplicaNumber( 0 );
    if ( theTouchable->GetHistoryDepth() > 0 ) {
      motherCopyNo = theTouchable->GetReplicaNumber( 1 );
    }
  }
  G4cout << "CCaloSD createNewHit for"
         << " PV "     << CurrentPV->GetName()
         << " PVid = " << currentCopyNo
         << " MVid = " << motherCopyNo
         << " Unit "   << UnitID <<G4endl;
  G4cout << " primary "    << PrimaryID
         << " time slice " << TSliceID 
         << " For Track  " << theTrack->GetTrackID()
         << " which is a " <<  theTrack->GetDefinition()->GetParticleName();
           
  if (theTrack->GetTrackID()==1) {
    G4cout << " of energy "     << theTrack->GetTotalEnergy();
  } else {
    G4cout << " daughter of part. " << theTrack->GetParentID();
  }

  G4cout  << " and created by " ;
  if (theTrack->GetCreatorProcess()!=NULL)
    G4cout << theTrack->GetCreatorProcess()->GetProcessName() ;
  else 
    G4cout << "NO process";
  G4cout << G4endl;
#endif          
    
  CurrentHit = new CCalG4Hit;
  CurrentHit->setTrackID(PrimaryID);
  CurrentHit->setTimeSlice(TSlice);
  CurrentHit->setUnitID(UnitID);
  CurrentHit->setEntry(EntrancePoint);
  CurrentHit->setIncidentEnergy(IncidentEnergy);
  updateHit();
  
  StoreHit(CurrentHit);

}         

void CCaloSD::updateHit() {
  if (EdepositEM+EdepositEHAD != 0) {
    CurrentHit->addEnergyDeposit(EdepositEM,EdepositEHAD);
#ifdef debug
    G4cout << "Energy deposit in Unit " << UnitID << " em " << EdepositEM/MeV
         << " hadronic " << EdepositEHAD/MeV << " MeV" << G4endl;
#endif
  }

  // buffer for next steps:
  TSID = TSliceID;
  PrimID = PrimaryID;
  PreviousPV = CurrentPV;
  PreviousUnitID = UnitID;
}

G4ThreeVector CCaloSD::SetToLocal(const G4ThreeVector& global) const{

  G4ThreeVector localPoint;
  const G4VTouchable*   touch= PreStepPoint->GetTouchable();
  localPoint=touch->GetHistory()->GetTopTransform().TransformPoint(global);
  
  return localPoint;  
}
     
void CCaloSD::EndOfEvent(G4HCofThisEvent*) {
  summarize();
}
     
void CCaloSD::summarize() {
}

void CCaloSD::clear() {
} 

void CCaloSD::DrawAll() {
} 

void CCaloSD::PrintAll() {
  G4cout << "CCaloSD: Collection " << theHC->GetName() << G4endl;
  theHC->PrintAllHits();
} 

void CCaloSD::SetOrganization(CCalVOrganization* org){

  if (theDescription!=0) 
    delete theDescription;
  theDescription = org;
}

G4double CCaloSD::curve_LY(const G4StepPoint* stepPoint) {

  G4double weight = 1.;
  G4ThreeVector localPoint = SetToLocal(stepPoint->GetPosition());
  const G4double crlength = 230.;
  G4double dapd = 0.5 * crlength - localPoint.z();
  if (dapd >= -0.1 || dapd <= crlength+0.1) {
    if (dapd <= 100.)
      weight = 1.05 - dapd * 0.0005;
  } else {
    G4cout << "CCaloSD, light coll curve : wrong distance to APD " << dapd 
           << " crlength = " << crlength
           << " z of localPoint = " << localPoint.z() 
           << " take weight = " << weight << G4endl;
  }
#ifdef ddebug
  G4cout << "CCaloSD, light coll curve : " << dapd 
         << " crlength = " << crlength
         << " z of localPoint = " << localPoint.z() 
         << " take weight = " << weight << G4endl;
#endif
  return weight;
}
