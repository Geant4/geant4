#ifndef Tst14PhysicsList_h
#define Tst14PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst14PhysicsListMessenger;
class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class Tst14PhysicsList: public G4VUserPhysicsList
{
public:

  Tst14PhysicsList();
  virtual ~Tst14PhysicsList();
  
protected:

  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
protected:

  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  
protected:
  // these methods Construct physics processes and register them
  void ConstructEM();
  void ConstructGeneral();

  void SetCuts();

public:

  void SetGammaCut(G4double);
  void SetElectronCut(G4double);

  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);
  void SetGELowLimit(G4double);
  void SetLowEnSecPhotCut(G4double);
  void SetLowEnSecElecCut(G4double);

private:

  G4LowEnergyIonisation*  LeIoprocess;
  G4LowEnergyPhotoElectric* LePeprocess;
  G4LowEnergyBremsstrahlung* LeBrprocess;
  Tst14PhysicsListMessenger* physicsListMessenger;
  G4double cutForGamma;
  G4double cutForElectron;

};

#endif







