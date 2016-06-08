//    ********************************************
//    *                                          *
//    *      BrachyPrimaryGeneratorAction.hh     *
//    *                                          *
//    ********************************************


#ifndef BrachyPrimaryGeneratorAction_h
#define BrachyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"

class G4ParticleGun;
class G4Event;

class BrachyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
	BrachyPrimaryGeneratorAction();
      ~BrachyPrimaryGeneratorAction();

 public:
      void GeneratePrimaries(G4Event* anEvent);

 private:
      G4ParticleGun* m_pParticleGun;
	G4RadioactiveDecay *m_pRadioactiveDecay;
};

#endif


