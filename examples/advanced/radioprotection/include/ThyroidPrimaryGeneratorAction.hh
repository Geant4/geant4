//
//    ********************************************
//    *                                          *
//    *      ThyroidPrimaryGeneratorAction.hh    *
//    *                                          *
//    ********************************************


#ifndef ThyroidPrimaryGeneratorAction_h
#define ThyroidPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
//class ThyroidAnalysisManager;
class ThyroidPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
	ThyroidPrimaryGeneratorAction();
      ~ThyroidPrimaryGeneratorAction();

 public:
      void GeneratePrimaries(G4Event* anEvent);

 private:
        G4ParticleGun* m_pParticleGun;
	G4RadioactiveDecay *m_pRadioactiveDecay;
        G4double Energy;     
};

#endif





