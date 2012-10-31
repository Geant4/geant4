#ifndef SimG4Core_PhysicsLists_QGSPCMS_FTFP_BERT_EML95msc93_H
#define SimG4Core_PhysicsLists_QGSPCMS_FTFP_BERT_EML95msc93_H

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
 
class QGSPCMS_FTFP_BERT_EML95msc93: public G4VModularPhysicsList {

public:

  QGSPCMS_FTFP_BERT_EML95msc93();
  virtual ~QGSPCMS_FTFP_BERT_EML95msc93();

  virtual void SetCuts();

};
 
#endif


