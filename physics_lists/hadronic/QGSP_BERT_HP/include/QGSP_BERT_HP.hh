#ifndef TQGSP_BERT_HP_h
#define TQGSP_BERT_HP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TQGSP_BERT_HP: public T
{
public:
  TQGSP_BERT_HP();
  virtual ~TQGSP_BERT_HP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "QGSP_BERT_HP.icc"
typedef TQGSP_BERT_HP<G4VModularPhysicsList> QGSP_BERT_HP;

// 2002 by J.P. Wellisch

#endif



