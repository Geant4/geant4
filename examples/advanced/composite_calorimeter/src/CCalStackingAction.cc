///////////////////////////////////////////////////////////////////////////////
// File: CCalStackingAction.cc
// Description: Stacking action needed for the application
///////////////////////////////////////////////////////////////////////////////
#include "CCalStackingAction.hh"
#include "G4StackManager.hh"

#include "G4SDManager.hh"
#include "CCaloSD.hh"
#include "SDList.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"

//#define debug
//#define ddebug

CCalStackingAction::CCalStackingAction():
isInitialized(false)
{}


CCalStackingAction::~CCalStackingAction(){}


void CCalStackingAction::PrepareNewEvent(){

  if(!isInitialized)   initialize();
  stage = firstStage;
  nurgent = 0;
  acceptSecondaries = 1;

}

void CCalStackingAction::initialize(){

  isInitialized = true;
 
  numberOfSD = SDList::getInstance()->getNumberOfCaloSD();
#ifdef debug
  cout << "CCalStackingAction look for " << numberOfSD 
       << " calorimeter-like SD" << endl;
#endif
  int i = 0;
  for (i=0; i<numberOfSD; i++) {
    G4String theName(SDList::getInstance()->getCaloSDName(i));
    SDName[i] = theName;
#ifdef debug
    cout << "Found SD  name " << theName << endl;
#endif
    theCaloSD[i] = 0;
  }   

  G4SDManager* sd = G4SDManager::GetSDMpointerIfExist();
  if (sd != 0) {
    
    for (i=0; i<numberOfSD; i++){

      G4VSensitiveDetector* aSD = sd->FindSensitiveDetector(SDName[i]);
      if (aSD==0) {
#ifdef debug
	cout << "CCalStackingAction::initialize: No SD with name " << SDName[i]
	     << " in this Setup " << endl;
#endif
      } else {
	theCaloSD[i] = dynamic_cast<CCaloSD*>(aSD);
	theCaloSD[i]->SetPrimaryID(0);
      }	   
    }
#ifdef debug
    cout << "CCalStackingAction::initialize: Could not get SD Manager !" 
	 << endl;
#endif
  }
   
}

G4ClassificationOfNewTrack CCalStackingAction::ClassifyNewTrack(const G4Track* aTrack){

  G4ClassificationOfNewTrack classification=fKill;
  int parentID = aTrack->GetParentID();
#ifdef ddebug
  G4TrackStatus status = aTrack->GetTrackStatus();
  cout << "Classifying track " << aTrack->GetTrackID()
       << " with status " << aTrack->GetTrackStatus() << endl;  
#endif
    
  if (parentID < 0) {
#ifdef debug     
    cout << "Killing track " << aTrack->GetTrackID() 
	 << " from previous event. Should not happen" << endl;
    cout << "returning classification= " << classification << endl;
#endif
    return classification= fKill;
  }
   
  if (aTrack->GetDefinition()->GetParticleName() == "gamma" && 
      aTrack->GetKineticEnergy() < 1.*eV) {
#ifdef debug
    cout << "Kills particle " << aTrack->GetDefinition()->GetParticleName() 
	 << " of energy " << aTrack->GetKineticEnergy()/MeV << " MeV" << endl;
#endif
    return classification= fKill;
  }
    
  if (stage<end) {
    /////////////////
    /// PRIMARIES ///
    /////////////////
    if (parentID == 0 ) {
      if ( nurgent == 0) {
	nurgent++;
	classification = fUrgent;
	setPrimaryID(aTrack->GetTrackID());
      }
      else  classification = fWaiting;   
    }

    ///////////////////
    /// SECONDARIES ///
    ///////////////////
       
    if (parentID > 0) {
      if (acceptSecondaries == 1) {
	if (trackStartsInCalo(const_cast<G4Track *>(aTrack))!=0 )
	  classification = fUrgent;
	else
	  classification = fWaiting; 
      } else {
	if(nurgent == 0){                     
	  nurgent++;
	  classification = fUrgent;
	  setPrimaryID(aTrack->GetTrackID());
	} else
	  classification = fWaiting;	
      }       
    }
       
       
  } else 
    classification = G4UserStackingAction::ClassifyNewTrack(aTrack);

#ifdef ddebug
  cout << " returning classification= " << classification
       << " for track "<< aTrack->GetTrackID() << endl;
#endif
  return classification;

}


void CCalStackingAction::NewStage(){

#ifdef ddebug
  cout << "In NewStage with stage = " << stage << endl;
#endif
  if (stage <end) {
    nurgent = 0;		    
    setPrimaryID(0);
    acceptSecondaries = 0;
    stackManager->ReClassify();
    acceptSecondaries = 1;
    if (stackManager->GetNUrgentTrack() == 0) {
      stage = stageLevel(stage+1);
    }
	
  }
}

G4bool CCalStackingAction::trackStartsInCalo(const G4Track* atrack){

 /// This method should check that the secondary particle
 /// was produced inside the detector calorimeter and 
 /// really is part of the shower.
 /// If it has been produced before the calorimeter 
 /// for ex. Bremsstrahlung, it should be treated as a new
 /// particle producing a new shower.

 return true;
}

void CCalStackingAction::setPrimaryID(G4int id){
  
  for (int i=0; i<numberOfSD; i++){
    if(theCaloSD[i] != 0)theCaloSD[i]->SetPrimaryID(id);
  }

}
