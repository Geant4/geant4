// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name: G4QMDNucleus.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 3 April 2007
// -------------------------------------------------------------------

#ifndef G4QMDNucleus_hh
#define G4QMDNucleus_hh

#include "G4QMDSystem.hh"
#include "G4QMDParameters.hh"

class G4QMDNucleus : public G4QMDSystem
{
   public:
      G4QMDNucleus();
      ~G4QMDNucleus();

      G4LorentzVector Get4Momentum();

      // Number of Nucleons (Proton or Neutron)
      G4int GetMassNumber();

      // Number of Protons
      G4int GetAtomicNumber();

      void CalEnergyAndAngularMomentumInCM();

      // rest mass from G4NucleiPropertiesTable
      G4double GetNuclearMass();

      void SetTotalPotential( G4double x ){ potentialEnergy = x; };
      G4double GetExcitationEnergy(){ return excitationEnergy; };

      G4int GetAngularMomentum(){ return jj; };

   private:

      G4double hbc;

      std::vector < G4ThreeVector > rcm, pcm;
      std::vector < G4double > es;
      G4int jj;

      G4double potentialEnergy;
      G4double excitationEnergy;
      G4double bindingEnergy;

      G4double kineticEnergyPerNucleon;
      //G4double bindingEnergyPerNucleon;
      //G4double potentialEnergyPerNucleon;
};

#endif
