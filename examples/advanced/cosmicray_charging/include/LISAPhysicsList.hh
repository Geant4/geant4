// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAPhysicsList_h
#define LISAPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class LISAPhysicsList: public G4VUserPhysicsList {

  public:
    LISAPhysicsList();
    virtual ~LISAPhysicsList();

  public:
    virtual void SetCuts();

  protected:
    // particles and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
    // physics processes
    virtual void AddTransportation();
    virtual void ElectromagneticPhysics();
    virtual void HadronicPhysics();
    virtual void ElectroNuclearPhysics();
    virtual void GeneralPhysics();

  private:
    G4int VerboseLevel;

};

#endif
