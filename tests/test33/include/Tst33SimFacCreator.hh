#ifndef Tst33SimFacCreator_hh
#define Tst33SimFacCreator_hh Tst33SimFacCreator_hh 

#include "Tst33VSimulationFactory.hh"


template<class S> class Tst33SimFacCreator : public Tst33VSimulationFactory {
public:
  Tst33SimFacCreator(const G4String &s)
    :
    fSimulationName(s)
  {}
  ~Tst33SimFacCreator(){}
  virtual const G4String &GetSimulationName() const{
    return fSimulationName;
  }
  virtual Tst33VSimulation *CreateSimulation() const{
    return new S;
  }
private:
  G4String fSimulationName;
};



#endif
