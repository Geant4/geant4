///////////////////////////////////////////////////////////////////////////////
// File: CCaloSD.cc
// Description: Stores hits of calorimetric type in appropriate container
///////////////////////////////////////////////////////////////////////////////

#include "CCaloSD.hh"
#include "G4VProcess.hh"
#include "G4SDManager.hh"
#include "G4VTouchable.hh"
#include "CCalVOrganization.hh"
#include "CCalSDList.hh"

#include "g4std/iostream"

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
 if (theDescription) 
   delete theDescription;
}


void CCaloSD::Initialize(G4HCofThisEvent*HCE) {

#ifdef debug
  G4cout << "CCaloSD : Initialize called for " << SDname << G4endl;
#endif
  //This initialization is performed at the beginning of an event
  //------------------------------------------------------------

  theHC = new CCalG4HitCollection(SDname, collectionName[0]);
  if (HCID<0) 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection( HCID, theHC );

  TSID = -2;
  PrimID = -2;
  /////PrimaryID = -99;  <--- initialized by StackingAction.
}


G4bool CCaloSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {

  if (aStep == NULL) return true;
   
  getStepInfo(aStep);
  if (hitExists() == false && EdepositEM+EdepositEHAD>0.) 
    createNewHit();

  return true;
} 


void CCaloSD::getStepInfo(G4Step* aStep) {
  
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


G4bool CCaloSD::hitExists() {
   
  if (PrimaryID<1) {
    G4cerr << "***** CCaloSD error: PrimaryID = " << PrimaryID
	   << " Maybe detector name changed"
	   << G4endl;
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

    CCalG4Hit* aPreviousHit = (*theHC)[j];
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
  G4cout << "CCaloSD createNewHit for"
       << " PV "     << CurrentPV->GetName()
       << " PVid = " << CurrentPV->GetCopyNo()
       << " MVid = " << CurrentPV->GetMother()->GetCopyNo()
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


G4ThreeVector CCaloSD::SetToLocal(G4ThreeVector global){

  G4ThreeVector localPoint;
  const G4VTouchable*   touch= PreStepPoint->GetTouchable();
  localPoint=touch->GetHistory()->GetTopTransform().TransformPoint(global);
  
  return localPoint;  

}
     

void CCaloSD::EndOfEvent(G4HCofThisEvent*HCE) {
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


