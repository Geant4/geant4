#ifndef TiaraSim_hh
#define TiaraSim_hh TiaraSim_hh
#include "G4UIterminal.hh"
#include "g4std/map"
#include <string>
#include "G4ProcessPlacer.hh"

class G4RunManager;
class G4VUserPrimaryGeneratorAction;
class G4VUserPhysicsList;
class TiaraPhysicsList;
class G4VUserDetectorConstruction;
class G4VisManager;
class TiaraStackingAction;
class G4UserEventAction;
class TiaraEnergyCutProcess;

class TiaraSim {
public:
  ~TiaraSim();
  static TiaraSim &GetTiaraSim();
  void SetGeometry(G4VUserDetectorConstruction *geo);
  void SetPrimaryGenerator(G4VUserPrimaryGeneratorAction *primGen);
  void SetPhysicsList(G4VUserPhysicsList *physList);
  void initialize();
  void initializeVisManager();
  void startSession();
  void AddParticleCut(const std::string &particle,
		      G4double cut);
  void AddTiaraEventAction(G4UserEventAction *eventAction);
  void AddVisRunAction();
  void BeamOn(G4int nEvents);

private:
  TiaraSim();
  static TiaraSim *fTiaraSim;
  G4RunManager *frunMgr;
  G4VisManager *fVisManager;
  G4UIterminal fSession;

  G4VUserDetectorConstruction *fGeometry;
  G4VUserPrimaryGeneratorAction *fPrimary;
  G4VUserPhysicsList *fPhysics;

  G4std::vector<TiaraEnergyCutProcess *> fCutProcessVector;
  G4std::vector<G4ProcessPlacer> fPlacers;

  TiaraStackingAction *fStackingAction;
  G4UserEventAction *fUserEventAction;
};




#endif










