// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AnalysisManager.cc,v 1.10 2000/11/21 15:41:25 barrand Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// Guy Barrand 20th Mai 2000
//
// Class description
//
//  This class is a plug-in for various analysis systems
// that follow the AIDA/Histogramming conventions.
//
#ifdef G4ANALYSIS_BUILD

#include "globals.hh"
#include "G4ios.hh"

#include "G4VAnalysisSystem.hh"
#include "G4AnalysisManager.hh"

G4AnalysisManager::G4AnalysisManager(){}
G4AnalysisManager::~G4AnalysisManager(){
  int number = fSystems.size();
  for(int count=0;count<number;count++) {
    delete fSystems[count];
  }
}
G4bool G4AnalysisManager::RegisterAnalysisSystem(
 G4VAnalysisSystem* aSystem
){
  fSystems.push_back(aSystem);
  return true;
}
G4VAnalysisSystem* G4AnalysisManager::FindSystem(const G4String& aName) {
  int number = fSystems.size();
  for(int count=0;count<number;count++) {
    if(aName==fSystems[count]->GetName()) return fSystems[count];
  }
  return 0;
}

IHistogramFactory* G4AnalysisManager::GetHistogramFactory(
 const G4String& aSystem
){
  G4VAnalysisSystem* analysisSystem = FindSystem(aSystem);
  if(!analysisSystem) {
    G4cerr << "G4AnalysisManager : " << G4std::endl;
    G4cerr << " " << aSystem << " unknown analysis system." << G4std::endl;
    return 0;
  }
  // Set it as current system.
  fCurrentSystem = analysisSystem;
  return analysisSystem->GetHistogramFactory();
}
void G4AnalysisManager::Store(IHistogram* aHistogram,const G4String& aSID){
  if(!fCurrentSystem) return;
  fCurrentSystem->Store(aHistogram,aSID);
}
void G4AnalysisManager::Plot(IHistogram* aHistogram){
  if(!fCurrentSystem) return;
  fCurrentSystem->Plot(aHistogram);
}

// Non AIDA things :
#ifdef G4ANALYSIS_NON_AIDA
ITuple* G4AnalysisManager::CreateTuple(
 const G4String& aSystem
,const G4String& aStorage
,const G4String& aSID
){
  G4VAnalysisSystem* analysisSystem = FindSystem(aSystem);
  if(!analysisSystem) {
    G4cerr << "G4AnalysisManager : " << G4std::endl;
    G4cerr << " " << aSystem << " unknown analysis system." << G4std::endl;
    return 0;
  }
  // Set it as current system.
  fCurrentSystem = analysisSystem;
  return analysisSystem->CreateTuple(aStorage,aSID);
}
ICloudFactory* G4AnalysisManager::GetCloudFactory(
 const G4String& aSystem
){
  G4VAnalysisSystem* analysisSystem = FindSystem(aSystem);
  if(!analysisSystem) {
    G4cerr << "G4AnalysisManager : " << G4std::endl;
    G4cerr << " " << aSystem << " unknown analysis system." << G4std::endl;
    return 0;
  }
  // Set it as current system.
  fCurrentSystem = analysisSystem;
  return analysisSystem->GetCloudFactory();
}
void G4AnalysisManager::Plot(ICloud* aCloud){
  if(!fCurrentSystem) return;
  fCurrentSystem->Plot(aCloud);
}
#endif

#endif
