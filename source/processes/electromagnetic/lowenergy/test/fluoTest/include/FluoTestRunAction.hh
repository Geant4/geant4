//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestRunAction_h
#define FluoTestRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "g4std/map"

#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisManager.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class FluoTestDataSet;
class G4DataVector;

class FluoTestRunAction : public G4UserRunAction
{
  public:
   
#ifdef G4ANALYSIS_USE
  FluoTestRunAction(FluoTestAnalysisManager* analysisMgr);
 FluoTestRunAction();
#else 
   FluoTestRunAction();
#endif 
 ~FluoTestRunAction();
  const FluoTestDataSet* GetSet();
  const FluoTestDataSet* GetGammaSet();
  const FluoTestDataSet* GetAlphaSet();
  const FluoTestDataSet* GetEfficiencySet();
  G4DataVector* GetEnergies();
  G4DataVector* GetData();
  

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  G4double GetDataSum();
  G4double GetInfData(G4double energy, G4double random);
  G4double GetSupData(G4double energy, G4double random);

private:
#ifdef G4ANALYSIS_USE
    FluoTestAnalysisManager* analysisManager;
#endif
  const FluoTestDataSet* dataSet;
  const FluoTestDataSet* dataGammaSet;
  const FluoTestDataSet* dataAlphaSet;
  const FluoTestDataSet* efficiencySet;
  G4DataVector* energies;
   G4DataVector* data;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > dataMap;
  
 void ReadData(G4double,G4String);
 void ReadResponse(const G4String& fileName);
  
 

};

#endif

