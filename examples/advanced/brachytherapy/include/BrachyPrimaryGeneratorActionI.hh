//    ********************************************
//    *                                          *
//    *      BrachyPrimaryGeneratorAction.hh     *
//    *                                          *
//    ********************************************


#ifndef BrachyPrimaryGeneratorActionI_h
#define BrachyPrimaryGeneratorActionI_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
#include "BrachyPrimaryGeneratorAction.hh"
class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyPrimaryGeneratorAction;
//class BrachyPrimaryGeneratorMessenger;
class BrachyPrimaryGeneratorActionI : public  G4VUserPrimaryGeneratorAction
{
 public:
	BrachyPrimaryGeneratorActionI();
      ~BrachyPrimaryGeneratorActionI();

 public:
      void GeneratePrimaries(G4Event* anEvent);
  
  G4double GetEnergy();
 private:
      G4ParticleGun* m_pParticleGun;
    
	G4RadioactiveDecay *m_pRadioactiveDecay;
        G4double Energy;
       
  //    BrachyPrimaryGeneratorMessenger* gunMessenger; 
   G4std::vector<G4double> vettore;
};

#endif


