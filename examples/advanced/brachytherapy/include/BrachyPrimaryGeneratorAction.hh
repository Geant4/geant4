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
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
	BrachyPrimaryGeneratorAction();
       ~BrachyPrimaryGeneratorAction();

 public:
  virtual  void GeneratePrimaries(G4Event* anEvent)=0;
  
     
 
    
};

#endif


