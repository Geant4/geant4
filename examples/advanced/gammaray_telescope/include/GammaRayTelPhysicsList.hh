//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelPhysicsList_h
#define GammaRayTelPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPhysicsList: public G4VUserPhysicsList
{
  public:
    GammaRayTelPhysicsList();
   ~GammaRayTelPhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
 
    virtual void SetCuts();

  public:
    // Set/Get cut values 
    void      SetCutForGamma(G4double);
    void      SetCutForElectron(G4double);
    void      SetCutForProton(G4double);           
    G4double  GetCutForGamma() const;
    G4double  GetCutForElectron() const;
    G4double  GetCutForProton() const;
    
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();

  private:
    G4double cutForGamma;
    G4double cutForElectron; 
    G4double cutForProton;
    G4double currentDefaultCut;
};

#endif



