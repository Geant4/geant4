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
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestRunAction::~hTestRunAction()
{
  delete hbookManager;
  histo.clear(); 
  G4cout << "runMessenger and histograms are deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start" << G4endl;
  zend     = 0.0;
  zend2    = 0.0;
  zEvt     = 0.0;
  verbose  = theDet->GetVerbose(); 
  nHisto   = theDet->GetHistoNumber(); 
  histName = theDet->GetHistoName();
  
#ifdef G4VIS_USE
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  }
#endif

  if(0 < nHisto) bookHisto();

  if(verbose > 0) {
    G4cout << "hTestRunAction: Histograms are booked and run has been started" 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::EndOfRunAction(const G4Run* aRun)
{

  G4cout << "hTestRunAction: End of run actions are started" << G4endl;

  // Zend average
  if(zEvt > 0.0) {
    zend  /= zEvt;
    zend2 /= zEvt;
    G4double sig = sqrt(zend2 - zend*zend);
    zend2 = sig / sqrt(zEvt);
    G4cout<<"========================================================="<<G4endl;
    G4cout << setprecision(4) << "Range(mm)= " << zend/mm 
           << "; Stragling(mm)= " << sig/mm 
           << setprecision(2) << " +- " << zend2/mm << G4endl;
    G4cout<<"========================================================="<<G4endl;
  }  

#ifdef G4VIS_USE
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
#endif

   // Write histogram file
  if(0 < nHisto) {
    hbookManager->write();
    G4cout << "Histograms and Ntuples are saved" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::SaveEvent()
{
  if(0 < nHisto) ntup->dumpData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::SaveToTuple(G4String parname, G4double val)
{
  if(0 < nHisto) ntup->column(parname,val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::SaveToTuple(G4String parname,G4double val,G4double defval)
{
  if(0 < nHisto) ntup->column(parname,val,defval);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::bookHisto()
{
  G4cout << "hTestRunAction: Histograms will be saved to the file <" 
         << histName << ">" << G4endl;

  // init hbook
  hbookManager = new HBookFile(histName, 68);

  // book histograms

  histo.resize(nHisto);

  G4int nbin = theDet->GetNumberOfAbsorbers();
  G4double zmax = (theDet->GetAbsorberThickness()) * G4double(nbin) / mm;

  if(0 < nHisto) histo[0] = hbookManager->histogram(
                "Energy deposit (MeV) in absorber (mm)",nbin,0.0,zmax);

  // book ntuple
  ntup = hbookManager->ntuple("Range/Energy");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::AddEnergy(G4double edep, G4double z)
{
  if(0 < nHisto) histo[0]->accumulate(z/mm, edep/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestRunAction::AddEndPoint(G4double z)
{
  zend  += z;
  zend2 += z*z;
  zEvt  += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....













