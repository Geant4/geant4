#ifndef MLPhysicsList_h
#define MLPhysicsList_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "MLVModularPhysicsList.hh"
#include "globals.hh"
#include <vector>

#ifndef USEHBOOK
  #include "RPTofstream.hh"
#endif

class MLPhysicsListMessenger;
class MLGeometryConstruction;
/*
class G4Regin;

struct Cuts {
  G4String regName;
  G4double gamma;
  G4double electron;
  G4double positron;
};
typedef  std::vector<Cuts> RegCuts;

*/
////////////////////////////////////////////////////////////////////////////////
//
class MLPhysicsList: public MLVModularPhysicsList
{
public:
  MLPhysicsList();
  virtual ~MLPhysicsList();
  
public:

  // Simulation scenario
  void SetScenario (G4String) ;
  G4String GetScenario () const {return scenario;};
  void ShowScenario () {
    G4cout <<"   Simulation Scenario: " <<scenario <<G4endl;};

  // Regions
  void AddARegion(G4String);
  void DeleteARegion(G4String);
  void ListRegions();

  void AddALayerToARegion(G4String, G4int);
  void DeleteALayerFromARegion(G4String, G4int);
  void ListTheLayersInARegion(G4String);
  
  // SetCuts() 
  virtual void SetCuts ();

  void SetGlobalDefault(G4double d) { cutForPositron 
					= cutForElectron = cutForGamma 
					=  defaultCutValue = d;};
  void SetGammaDefault (G4double d) {cutForGamma = d;};
  void SetElectronDefault (G4double d) {cutForElectron = d;};
  void SetPositronDefault (G4double d) {cutForPositron = d;};

  void SetCutInRangeForRegion (G4double, G4String);
  void SetCutInRangeForRegionForParticle( G4double, G4String, G4String);


private:

  void BuildList ();

  void SetGlobalCuts ();
  void SetCutsByRegion ();

private:

  MLPhysicsListMessenger* physListMessenger;
  
  G4String scenario;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  //  std::vector<G4Regin*> mlRegion;
  //  RegCuts  mlRegCuts;

  MLGeometryConstruction* geometry;

#ifndef USEHBOOK
  friend RPTofstream & operator << (RPTofstream &, const MLPhysicsList &);
#endif
};
////////////////////////////////////////////////////////////////////////////////
#endif
