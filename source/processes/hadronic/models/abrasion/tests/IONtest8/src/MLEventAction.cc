////////////////////////////////////////////////////////////////////////////////
//
#include "MLEventAction.hh"
#include "MLHit.hh"
#include "MLEventActionMessenger.hh"
#include "MLRunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4DynamicParticle.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"


// This file is a global variable in which we store energy deposition per hit
// and other relevant information
// extern std::ofstream outFile;
////////////////////////////////////////////////////////////////////////////////
//
MLEventAction::MLEventAction(MLGeometryConstruction* det)
  :edetectorID(-1), geometry(det), drawFlag("all"),printModulo(100),cputime(600.)
{
  eventMessenger  = new MLEventActionMessenger(this);
  analysisManager = MLAnalysisManager::getInstance();
  fluenceAnalyser = analysisManager->GetFluenceAnalyser();
  NIELAnalyser    = analysisManager->GetNIELAnalyser();
  doseAnalyser    = analysisManager->GetDoseAnalyser();
  PHSAnalyser     = analysisManager->GetPHSAnalyser();
}
////////////////////////////////////////////////////////////////////////////////
//
MLEventAction::~MLEventAction()
{
  delete eventMessenger;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLEventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  if (evtNb == 0) {
    
    G4cout <<G4endl <<"Total CPU time used for initialization: "
           <<(static_cast <MLRunManager*> (MLRunManager::GetRunManager()))
             ->GetTimeUsed()
           <<" seconds" <<G4endl;
  }
  if (!((evtNb+1)%printModulo)) {
    G4cout <<G4endl <<"Event: "<<evtNb+1 <<G4endl;
    HepRandom::showEngineStatus();
  }
  if (edetectorID == -1) {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    edetectorID = SDman->GetCollectionID("MLCollection");
  }
  NIELAnalyser->BeginOfEventAction();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id       = evt->GetEventID();
  G4int NbOfLayers     = geometry->GetNbOfLayers();
  G4int NbOfELayers    = geometry->GetNbOfELayers();
  G4int NbOfDLayers    = doseAnalyser->GetNbOfDLayers();
  G4int NbOfNielLayers = NIELAnalyser->GetNbOfNielLayers();
  
  G4HCofThisEvent* HCE  = evt->GetHCofThisEvent();
  MLHitsCollection* EHC = NULL;
  
  if (HCE) EHC = (MLHitsCollection*)(HCE->GetHC(edetectorID));
  
  if (EHC) {
    G4int n_hit        = EHC->entries();
    G4double* Edeposit = new G4double[NbOfLayers];
    G4double* Eweight  = new G4double[NbOfLayers];
    G4int i, j;
    G4int k            = 0;
    for ( j=0; j<NbOfLayers; j++) {Edeposit[j] = Eweight[j] = 0.;};

    for ( i=0; i<n_hit; i++) {
      G4int layer = (*EHC)[i]->GetLayer();
      G4double ener  = (*EHC)[i]->GetEdep();
      G4double wei   = (*EHC)[i]->GetWeight();
      if (ener > 0.) {
	//        layeridx            = geometry->GetELayerPosition(layer);
        Edeposit[layer] += ener;
	Eweight[layer] += (wei*ener);
      } else {
        G4ThreeVector aparticle = (*EHC)[i]->GetParticle();
        G4double weight         = (*EHC)[i]->GetWeight();
        G4int pName             = (G4int) aparticle.x();
        G4double energy         = aparticle.y();
        G4double theta          = aparticle.z();
	//	  layeridx = geometry->GetFLayerPosition(layer);
	//cout <<  pName << " " << energy << " " << theta << " " << layer <<  endl;
        fluenceAnalyser->TallyFluenceEvent(layer,pName,theta,energy,weight);
        for (j=0; j<NbOfNielLayers; j++) {
          k = NIELAnalyser->GetNLayerIdx(j);
          if (k == layer) {
            NIELAnalyser->AddToNIELDetector(j,pName,theta,energy,weight);
            break;
          }
        }
      }
    }
    for (j=0; j<NbOfELayers; j++) {
      // Here we fill the histograms of the Analysis manager
      k = geometry->GetELayerIdx(j);
      if (Edeposit[k] > 0.)  PHSAnalyser->FillHbook(j,Edeposit[k],0.,Eweight[k]/Edeposit[k]);
    }
    for (j=0; j<NbOfDLayers; j++) {
      k = doseAnalyser->GetDLayerIdx(j);
      if (Edeposit[k] > 0.) doseAnalyser->AddToDoseDetector(j,Eweight[k]);
    }

    delete [] Edeposit;
    delete [] Eweight;
  }

  NIELAnalyser->EndOfEventAction();
  // Here we call the analysis manager for persistency
  if ((event_id)%printModulo==0  && event_id != 0) {
    analysisManager->EndOfEventAction(G4double(event_id));
    if ((static_cast <MLRunManager*>
	 (MLRunManager::GetRunManager()))->GetTimeUsed() > cputime) {
      G4Exception(" Allocated CPU time exceeded." );
    } else {
//      HepRandom::showEngineStatus();
      G4cout <<"Total CPU time used so far: "
             <<(static_cast <MLRunManager*> (MLRunManager::GetRunManager()))
	->GetTimeUsed()
             <<" seconds." <<G4endl;
    }
  }
  if (event_id < 100 && G4VVisManager::GetConcreteInstance()) {
    G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    for (G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory *)
        ((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") {
        trj->DrawTrajectory(0);
      } else if (drawFlag == "charged" && trj->GetCharge() != 0.) {
        trj->DrawTrajectory(0);
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLEventAction::SetRandomSeed(G4int index = 0)
{
  const  G4long* seeds = HepRandom::getTheSeeds ();
  HepRandom::setTheSeeds(seeds,index);
}
////////////////////////////////////////////////////////////////////////////////
