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
//      ---------- hTestCalorimeterSD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestRunAction.hh"
#include "hTestRunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

#include "CLHEP/Hist/HBookFile.h"
#include <assert.h>

#include "hTestPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestRunAction::hTestRunAction(hTestDetectorConstruction* det):
  theDet(det),
  histName("histo.hbook"),
  hbookManager(0),
  ntup(0),
  nHisto(1),
  verbose(0)
{
  runMessenger = new hTestRunMessenger(this);

  // Water test
  // nbinEn = 600;
  // Enlow  = 0.0;
  //Enhigh = 60.0;
  // Mylar/GaAs test
  //  nbinEn = 60;
  //Enlow  = 0.0;
  //Enhigh = 0.06;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestRunAction::~hTestRunAction()
{
  delete runMessenger;
  delete hbookManager;
  histo.clear(); 
  G4cout << "runMessenger and histograms are deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  verbose = theDet->GetVerbose(); 
  histName = theDet->GetHistoName();
  G4cout << "### Run " << aRun->GetRunID() << " start" << G4endl;
  
#ifdef G4VIS_USE
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  }
#endif

  bookHisto();

  G4cout << "hTestRunAction: Histograms are booked" << G4endl;
      
  G4cout << "hTestRunAction: Run is started!!!" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::EndOfRunAction(const G4Run* aRun)
{

  G4cout << "hTestRunAction: End of run actions are started" << G4endl;

#ifdef G4VIS_USE
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
#endif

   // Write histogram file
  hbookManager->write();
  G4cout << "Histograms and Ntuples are saved" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::SaveEvent()
{
  ntup->dumpData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::SaveToTuple(G4String parname, G4double val)
{
  ntup->column(parname,val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::SaveToTuple(G4String parname,G4double val,G4double defval)
{
  ntup->column(parname,val,defval);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::bookHisto()
{
  G4cout << "hTestRunAction: Histograms will be saved to the file <" 
         << histName << ">" << G4endl;

  // init hbook
  hbookManager = new HBookFile(histName, 68);

  // book histograms

  histo.resize(nHisto)

  G4int nbin = theDet->GetAbsorberNumber();
  G4int zmax = (theDet->GetAbsorberThickness()) * nbin / mm;

  if(0 < nHisto) {
    histo[0] = hbookManager->histogram("Energy deposit (MeV) in absorber (mm)"
                                       ,nbin,0.0,zmax);
  }

  // book ntuple
  ntup = hbookManager->ntuple("Range/Energy");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::AddEnergy(G4double edep, G4double z)
{
  if(0 < nHisto) {
    histo[0]->accumulate(z/mm, edep/MeV);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....













