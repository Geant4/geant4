// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LabSystem.cc,v 1.2 2000/11/16 13:44:41 barrand Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 14 September 2000

#ifdef G4ANALYSIS_BUILD_LAB

#include <Lab/Analysis.h>

#include "G4ios.hh"

#include "G4LabSystem.hh"

#ifdef G4ANALYSIS_LAB_VISUALIZATION
#include <ISession.h>
#include <IUI.h>
#include <IViewer.h>
#include <Lib/Compiler.h>
BegNameSpace(Lib)
#include <Lib/Memory.h>
EndNameSpace
BegNameSpace(OnX)
#include <OnX/Interpreters.h>
#include <OnX/Session.h>
#include <OnX/XmUI.h>
EndNameSpace
#endif

G4LabSystem::G4LabSystem (
 const G4String& aName
)
:fName(aName)
,fStorageManager(0)
,fHistogramManager(0)
,fCloudManager(0)
,fTupleManager(0)
,fSession(0)
{
#ifdef G4ANALYSIS_LAB_VISUALIZATION
  fSession = new OnX::Session;
  if(fSession) {
    OnX::XmUI* gui = new OnX::XmUI(fSession,0,0);
    if(gui) {
      fSession->addManager(gui);
      gui->executeCreateCallbacks();
      //gui->steer();
    }
  }
  fHistogramManager = 
    (IHistogramManager*)fSession->findManager("HistogramManager");
  if(!fHistogramManager) {
    G4cout << "Can't find the histogram manager" << G4std::endl;
  }
  fCloudManager = 
    (ICloudManager*)fSession->findManager("CloudManager");
  if(!fCloudManager) {
    G4cout << "Can't find the cloud manager" << G4std::endl;
  }
  fTupleManager = 
    (ITupleManager*)fSession->findManager("TupleManager");
  if(!fTupleManager) {
    G4cout << "Can't find the tuple manager" << G4std::endl;
  }
  fStorageManager = 
    (IStorageSystemManager*)fSession->findManager("StorageManager");
  if(!fStorageManager) {
    G4cout << "Can't find the storage manager" << G4std::endl;
  }
#else
  fHistogramManager = new Lab::HistogramManager();
  fCloudManager = new Lab::CloudManager();
  fTupleManager = new Lab::TupleManager();
  fStorageManager = new Lab::StorageManager();
#endif
}
G4LabSystem::~G4LabSystem (){
#ifdef G4ANALYSIS_LAB_VISUALIZATION
  delete fSession; // Will delete the gui object and the managers.
#else
  delete fHistogramManager;
  delete fCloudManager;
  delete fTupleManager;
  delete fStorageManager;
#endif
  Lib::Memory::dump();
}
const G4String& G4LabSystem::GetName() const {
  return fName;
}
IHistogramFactory* G4LabSystem::GetHistogramFactory() {
  return fHistogramManager;
}
ICloudFactory* G4LabSystem::GetCloudFactory() {
  return fCloudManager;
}
ITuple* G4LabSystem::CreateTuple(
 const G4String& aStorage
,const G4String& aSID
){
  if(!fTupleManager) return 0;
  Lab::RioStorage* storageRio = new Lab::RioStorage(aStorage,
						    "RECREATE",
						    fStorageManager);
  return new Lab::RioTuple(storageRio,aSID,fTupleManager);
}

void G4LabSystem::Store(IHistogram*,const G4String&){
  if(!fHistogramManager) return;
  Lab::RioStorage* storage = new Lab::RioStorage("g4osc.root","RECREATE");
  //rioDB->createDirectory("histograms");
  //rioDB->cd("histograms");
  int number = fHistogramManager->size();
  for(int count=0;count<number;count++) {
    IHistogram* histo = (*fHistogramManager)[count];
    const char* title = histo->title().c_str();
    //G4cout << "G4Lab : store : " << title << G4std::endl;
    fHistogramManager->store(histo,storage,title);
  }
  storage->commit();
  delete storage;
  //if(verboseLevel>1) G4cout << "Analysis is deleting." << G4endl;
  //clear();
}
void G4LabSystem::Plot(IHistogram* aHistogram){
  if(!fHistogramManager) return;
#ifdef G4ANALYSIS_LAB_VISUALIZATION
  if(fSession) {
    IUI* ui = (IUI*)fSession->findManager("UI_Manager");
    if(ui) {
      IViewer* viewer = ui->getCurrentViewer();
      if(viewer) {
	viewer->erase("region");
	fHistogramManager->visualize(aHistogram,viewer);
	viewer->next();
	ui->synchronize();
      }
    }
  }
#endif
}
void G4LabSystem::Plot(ICloud* aCloud){
  if(!fCloudManager) return;
#ifdef G4ANALYSIS_LAB_VISUALIZATION
  if(fSession) {
    IUI* ui = (IUI*)fSession->findManager("UI_Manager");
    if(ui) {
      IViewer* viewer = ui->getCurrentViewer();
      if(viewer) {
	viewer->erase("region");
	fCloudManager->visualize(aCloud,viewer);
	viewer->next();
	ui->synchronize();
      }
    }
  }
#endif
}

#endif
