#ifdef G4ANALYSIS_USE
#ifndef XrayFluoAnalysisManager_h
#define XrayFluoAnalysisManager_h 1

#include "globals.hh"

// Histogramming from AIDA 
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// Histogramming from Anaphe
#include "Interfaces/IHistoManager.h"

// Ntuples from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
#include "NtupleTag/LizardQuantity.h"


class G4Step;
class XrayFluoAnalysisMessenger;
class NTuple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoAnalysisManager
{
public:
 
  virtual ~XrayFluoAnalysisManager();
  
  void book();
  
  void finish();
  
  void analyseStepping(const G4Step* aStep);
 
  void analyseEnergyDep(G4double eDep);

 void analysePrimaryGenerator(G4double energy);

 static XrayFluoAnalysisManager* getInstance();
 

private:
  
  XrayFluoAnalysisManager();
 
  static XrayFluoAnalysisManager* instance;
  
  IHistoManager* histoManager;
 
  Lizard::NTupleFactory* factory;
  Lizard::NTuple* ntuple;
  
  XrayFluoAnalysisMessenger* analysisMessenger;

  // Quantities for the ntuple
  Lizard::Quantity<float> eDep;
  Lizard::Quantity<float> counts;
};
#endif
#endif



