// -------------------------------------------------------------------
// $Id: MicrobeamPhysicsList.hh,v 1.1 2006-04-06 15:32:43 sincerti Exp $
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



