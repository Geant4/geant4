
#ifndef exGPSPhysicsList_h
#define exGPSPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSPhysicsList: public G4VUserPhysicsList
{
  public:
    exGPSPhysicsList();
   ~exGPSPhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  public:
     
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();
    void ConstructNuclei();

  private:
};

#endif
