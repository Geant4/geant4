#ifndef G4HadronicInteractionRegistry_h
#define G4HadronicInteractionRegistry_h 1

#include <rw/tpvector.h>
#include "globals.hh"
class G4HadronicInteraction;

class G4HadronicInteractionRegistry
{
  public:
  
  ~G4HadronicInteractionRegistry();
  G4HadronicInteractionRegistry(){nModels = 0;}
  
  static void RegisterMe(G4HadronicInteraction * aModel);
  static void RemoveMe(G4HadronicInteraction * aModel){};
  
  private:
  //  !!!  can not use "copy constructor" nor "default constructor" !!!!
       G4HadronicInteractionRegistry(const G4HadronicInteractionRegistry &right) 
       { nModels = right.nModels; }
       G4HadronicInteractionRegistry() {nModels = 0;}

  //  !!!  Assignment operation is forbidden !!!
      const G4HadronicInteractionRegistry & operator=(const G4HadronicInteractionRegistry &right) 
      { return *this;}

  void AddModel(G4HadronicInteraction * aModel);
  
  G4int nModels;
  RWTPtrVector<G4HadronicInteraction> allModels;
  static G4HadronicInteractionRegistry theRegistry;

};

#endif;
