// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestEventAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestEventAction.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::hTestEventAction(hTestRunAction* run, 
                                   hTestDetectorConstruction* det):
  theRun(run),
  theDet(det),
  verbose(0),
  nEvt(0),
  drawFlag("all")
{
  theDet->SetEventAction(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::~hTestEventAction()
{
  energy.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::BeginOfEventAction(const G4Event* evt)
{  
  verbose = theRun->GetVerbose();
  nEvt++;
  if(verbose > 0) {
    G4cout << "hTestEventAction: Event #" 
           << evt->GetEventID() << " started; nEvt = " 
           << nEvt << G4endl;
  }

  //theRun->InitializeTuples();

  numAbs = theDet->GetNumberOfAbsorbers();
  energy.resize(numAbs);
  for(G4int i=0; i<numAbs; i++) { energy[i] = 0.0; }

  backEnergy = 0.0;
  leakEnergy = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int i, j;

  theRun->SaveToTuple(G4String("backE"),backEnergy);      
  theRun->SaveToTuple(G4String("leakE"),leakEnergy);      

  // The histogramm on the energy deposition profile
  if(numAbs > 0) {
    G4double s = theDet->GetAbsorberThickness();
    G4double z = -0.5 * s;
    for(i=0; i<numAbs; i++) {
      z += s; 
      theRun->AddEnergy(energy[i], z);
    }
  }

  // Integrated energy deposition to nTuple
  G4int nMax = 60;
  G4double EE[60];
  G4String eSlice[60]={
      "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", 
      "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", 
      "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", 
      "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", 
      "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", 
      "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59"};
  G4String eInteg[60]={
      "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", 
      "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", 
      "E20", "E21", "E22", "E23", "E24", "E25", "E26", "E27", "E28", "E29", 
      "E30", "E31", "E32", "E33", "E34", "E35", "E36", "E37", "E38", "E39", 
      "E40", "E41", "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49", 
      "E50", "E51", "E52", "E53", "E54", "E55", "E56", "E57", "E58", "E59"};

  if (nMax > numAbs) nMax = numAbs;

  if(nMax > 1) {
    for(i=0; i<nMax; i++){
      EE[i]=0.0;
      for(j=0; j<i+1; j++) {
        EE[i] += energy[j];
      }  
      theRun->SaveToTuple(eSlice[i],energy[i]);      
      theRun->SaveToTuple(eInteg[i],EE[i]);      
    }
  }

#ifdef G4VIS_USE  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    G4TrajectoryContainer* trjc = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trjc) n_trajectories = trjc->entries();  

    for(i=0; i<n_trajectories; i++) {
      G4Trajectory* t = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") t->DrawTrajectory(50);
      else if ((drawFlag == "charged")&&(t->GetCharge() != 0.))
                             t->DrawTrajectory(50); 
    }
  }  
#endif

  if(verbose > 0) {
    G4cout << "hTestEventAction: Event #" 
           << evt->GetEventID() << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::AddEnergy(G4double edep, G4int n)
{  
  if(n < 0 || n >= numAbs) {
    G4cout << "Warning!!! hTestEventAction: cannot add " << edep/MeV
           << " MeV to the slice # " << n << G4endl;
    return;
  }

  energy[n] += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  













