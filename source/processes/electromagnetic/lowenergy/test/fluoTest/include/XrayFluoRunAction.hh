//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoRunAction_h
#define XrayFluoRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "g4std/map"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class XrayFluoDataSet;
class G4DataVector;

class XrayFluoRunAction : public G4UserRunAction
{
  public:
   

 XrayFluoRunAction();

 ~XrayFluoRunAction();
  const XrayFluoDataSet* GetSet();
  const XrayFluoDataSet* GetGammaSet();
  const XrayFluoDataSet* GetAlphaSet();
  const XrayFluoDataSet* GetEfficiencySet();
  G4DataVector* GetEnergies();
  G4DataVector* GetData();
  

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  G4double GetDataSum();
  G4double GetInfData(G4double energy, G4double random);
  G4double GetSupData(G4double energy, G4double random);

private:

  const XrayFluoDataSet* dataSet;
  const XrayFluoDataSet* dataGammaSet;
  const XrayFluoDataSet* dataAlphaSet;
  const XrayFluoDataSet* efficiencySet;
  G4DataVector* energies;
   G4DataVector* data;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > dataMap;
  
 void ReadData(G4double,G4String);
 void ReadResponse(const G4String& fileName);
  
 

};

#endif

