// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAPrimaryGeneratorAction_h
#define LISAPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class LISAPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

   public:
      LISAPrimaryGeneratorAction();
      ~LISAPrimaryGeneratorAction();
      void GeneratePrimaries(G4Event* anEvent);

   private:
      G4GeneralParticleSource* particleGun;

   private:
      long seeds[2];
      G4double energy_pri;

   public:
      const long* GetEventSeeds()   {return seeds;};
      G4double GetEnergyPrimary()   {return energy_pri;};

};

#endif

