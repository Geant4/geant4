#ifndef BrachyPrimaryGeneratorActionIr_h
#define BrachyPrimaryGeneratorActionIr_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
#include "BrachyPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyPrimaryGeneratorAction;
class BrachyPrimaryGeneratorActionIr: public G4VUserPrimaryGeneratorAction
{
 public:
	BrachyPrimaryGeneratorActionIr();
      ~BrachyPrimaryGeneratorActionIr();

 public:
      void GeneratePrimaries(G4Event* anEvent);
  G4double GetEnergy(){return Energy;};
     
 private:
  // G4ParticleGun* m_pParticleGun;
     G4ParticleGun* m_pParticleGun;
  //G4RadioactiveDecay *m_pRadioactiveDecay;
        G4double Energy;
       G4RadioactiveDecay *m_pRadioactiveDecay;
  //    BrachyPrimaryGeneratorMessenger* gunMessenger; 
  G4std::vector<G4double> vettore;
};

#endif






























































