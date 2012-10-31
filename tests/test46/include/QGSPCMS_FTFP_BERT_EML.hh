#ifndef SimG4Core_PhysicsLists_QGSPCMS_FTFP_BERT_EML_H
#define SimG4Core_PhysicsLists_QGSPCMS_FTFP_BERT_EML_H
 
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
 
class QGSPCMS_FTFP_BERT_EML: public G4VModularPhysicsList {

public:

  QGSPCMS_FTFP_BERT_EML();
  virtual ~QGSPCMS_FTFP_BERT_EML();

  virtual void SetCuts();

};

#endif


