//    **********************************
//    *                                *
//    *      BrachyPhysicsList.hh      *
//    *                                *
//    **********************************


#ifndef BrachyPhysicsList_h
#define BrachyPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class BrachyPhysicsList: public G4VUserPhysicsList
{
 public:
    BrachyPhysicsList();
    ~BrachyPhysicsList();

 protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();
  
 private:
    G4ParticleDefinition *m_pGamma;
    G4ParticleDefinition *m_pPositron;
    G4ParticleDefinition *m_pElectron;
 
   void ConstructEM();
};
#endif
