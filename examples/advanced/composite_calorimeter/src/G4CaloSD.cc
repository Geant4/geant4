///////////////////////////////////////////////////////////////////////////////
// File: G4CaloSD.cc
// Author: Veronique Lefebure
// Modifications: 27/04/00 S.B.
///////////////////////////////////////////////////////////////////////////////

#include "G4CaloSD.hh"
#include "G4VProcess.hh"
#include "G4SDManager.hh"
#include "G4VTouchable.hh"
#include "VDetectorOrganization.hh"
#include "SDList.hh"

#include<iostream>

//#define debug
//#define ddebug
 
G4CaloSD::G4CaloSD(G4String name, VDetectorOrganization* numberingScheme):
  G4VSensitiveDetector(name), SDname(name), theDescription(numberingScheme),
  CurrentHit(0), theTrack(0), CurrentPV(0), PreviousPV(0), PreStepPoint(0), 
  PostStepPoint(0), theHC(0), HCID(-1), UnitID(0), PreviousUnitID(0) {
  
  collectionName.insert(name);
  
  cout << "*******************************************************" << endl;
  cout << "*                                                     *" << endl;
  cout << "* Constructing a G4CaloSD  with name " << name           << endl;  
  cout << "*                                                     *" << endl;
  cout << "*******************************************************" << endl;

  SDList::getInstance()->addCalo(name);
}


G4CaloSD::~G4CaloSD() {
 if (theDescription) 
   delete[] theDescription;
}


void G4CaloSD::Initialize(G4HCofThisEvent*HCE) {

#ifdef debug
  cout << "G4CaloSD : Initialize called for " << SDname << endl;
#endif
  //This initialization is performed at the beginning of an event
  //------------------------------------------------------------

  theHC = new G4CaloHitsCollection(SDname, collectionName[0]);
  if (HCID<0) 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection( HCID, theHC );

  TSID = -2;
  PrimID = -2;
  /////PrimaryID = -99;  <--- initialized by StackingAction.
}


G4bool G4CaloSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {

  if (aStep == NULL) return true;
   
  getStepInfo(aStep);
  if (hitExists() == false && EdepositEM+EdepositEHAD>0.) 
    createNewHit();

  return true;
} 


void G4CaloSD::getStepInfo(G4Step* aStep) {
  
  PreStepPoint = aStep->GetPreStepPoint(); 
  PostStepPoint= aStep->GetPostStepPoint(); 
  HitPoint     = PreStepPoint->GetPosition();	

  theTrack = aStep->GetTrack();   
  CurrentPV= PreStepPoint->GetPhysicalVolume();

  G4String     particleType =  theTrack->GetDefinition()->GetParticleName();

  if (particleType == "e-" ||
      particleType == "e+" ||
      particleType == "gamma" ){
    EdepositEM   = aStep->GetTotalEnergyDeposit();
    EdepositEHAD = 0.;
  } else {
    EdepositEM   = 0.;
    EdepositEHAD = aStep->GetTotalEnergyDeposit();
  }

  TSlice = (PostStepPoint->GetGlobalTime() )/nanosecond;
  TSliceID = (int) TSlice;
  if (theDescription!=0) 
    UnitID = theDescription->GetUnitID(aStep);
  else
    UnitID = 0;
   
}


G4bool G4CaloSD::hitExists() {
   
  if (PrimaryID<1) {
    G4cerr << "***** G4CaloSD error: PrimaryID = " << PrimaryID
	   << " Maybe detector name changed"
	   << endl;
  }
   
      
  if ( CurrentPV==PreviousPV && PrimaryID == PrimID && TSliceID == TSID &&
       UnitID==PreviousUnitID) {
    updateHit();
    return true;
  }
   

  if (PrimaryID != PrimID)
    ResetForNewPrimary();
   

  //look in HC whether a hit with the same primID,UnitID,TSliceID exists already:
   
  G4bool found = false;
  for (int j=0; j<theHC->entries()&&!found; j++) {

    G4CaloHit* aPreviousHit = (*theHC)[j];
    if (aPreviousHit->getTrackID()  == PrimaryID &&
	aPreviousHit->getTimeSliceID() == TSliceID  &&
	aPreviousHit->getUnitID()== UnitID       ) {
      CurrentHit = aPreviousHit;
      found = true;
    }
  }          

  if (found) {
    updateHit();
    return true;
  } else {
    return false;
  }    
}


void G4CaloSD::ResetForNewPrimary() {
  
  EntrancePoint = SetToLocal(HitPoint);
  IncidentEnergy = PreStepPoint->GetKineticEnergy();

}


void G4CaloSD::StoreHit(G4CaloHit* hit){

  if (PrimID<0) return;
  if (hit == 0 ) {
    cout << "G4CaloSD: hit to be stored is NULL !!" <<endl;
    return;
  }

  theHC->insert( hit );
}


void G4CaloSD::createNewHit() {

#ifdef debug
  cout << "G4CaloSD createNewHit for"
       << " PV "     << CurrentPV->GetName()
       << " PVid = " << CurrentPV->GetCopyNo()
       << " MVid = " << CurrentPV->GetMother()->GetCopyNo()
       << " Unit "   << UnitID <<endl;
  cout << " primary "    << PrimaryID
       << " time slice " << TSliceID 
       << " For Track  " << theTrack->GetTrackID()
       << " which is a " <<  theTrack->GetDefinition()->GetParticleName();
	   
  if (theTrack->GetTrackID()==1) {
    cout << " of energy "     << theTrack->GetTotalEnergy();
  } else {
    cout << " daughter of part. " << theTrack->GetParentID();
  }

  cout  << " and created by " ;
  if (theTrack->GetCreatorProcess()!=NULL)
    cout << theTrack->GetCreatorProcess()->GetProcessName() ;
  else 
    cout << "NO process";
  cout << endl;
#endif          
    

  CurrentHit = new G4CaloHit;
  CurrentHit->setTrackID(PrimaryID);
  CurrentHit->setTimeSlice(TSlice);
  CurrentHit->setUnitID(UnitID);
  CurrentHit->setEntry(EntrancePoint);
  CurrentHit->setIncidentEnergy(IncidentEnergy);
  updateHit();
  
  StoreHit(CurrentHit);

}	 


void G4CaloSD::updateHit() {
  if (EdepositEM+EdepositEHAD != 0) {
    CurrentHit->addEnergyDeposit(EdepositEM,EdepositEHAD);
#ifdef debug
    cout << "Energy deposit in Unit " << UnitID << " em " << EdepositEM/MeV
	 << " hadronic " << EdepositEHAD/MeV << " MeV" << endl;
#endif
  }

  // buffer for next steps:
  TSID = TSliceID;
  PrimID = PrimaryID;
  PreviousPV = CurrentPV;
  PreviousUnitID = UnitID;
}


G4ThreeVector G4CaloSD::SetToLocal(G4ThreeVector global){

  G4ThreeVector localPoint;
  const G4VTouchable*   touch= PreStepPoint->GetTouchable();
  localPoint=touch->GetHistory()->GetTopTransform().TransformPoint(global);
  
  return localPoint;  

}
     

void G4CaloSD::EndOfEvent(G4HCofThisEvent*HCE) {
  summarize();
}
     

void G4CaloSD::summarize() {
}


void G4CaloSD::clear() {
} 


void G4CaloSD::DrawAll() {
} 


void G4CaloSD::PrintAll() {
  cout << "G4CaloSD: Collection " << theHC->GetName() << endl;
  theHC->PrintAllHits();
} 


void G4CaloSD::SetOrganization(VDetectorOrganization* org){

  if (theDescription!=0) 
    delete theDescription;
  theDescription = org;

}


