// -------------------------------------------------------------------
// $Id: MicrobeamPhysicsList.hh,v 1.2 2006-04-10 14:47:31 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamPhysicsList_h
#define MicrobeamPhysicsList_h 1

#include "G4VUserPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrobeamPhysicsList: public G4VUserPhysicsList
{
public:
  MicrobeamPhysicsList();
  ~MicrobeamPhysicsList();

private:
  
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  
protected:
  
  void ConstructParticle();
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructBaryons();

  void ConstructProcess();
  void ConstructEM();
  void ConstructHad();
  void ConstructGeneral();
  void SetCuts();
  
};
#endif



