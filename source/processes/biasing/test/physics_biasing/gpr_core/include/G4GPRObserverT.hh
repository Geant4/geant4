#ifndef G4GPROBSERVERT_HH
#define G4GPROBSERVERT_HH

#include"G4GPRObserverCollectionT.hh"
#include "G4GPRTriggerTypes.hh"

template <typename Type>
struct G4GPRObserverT : public Type::ObserverCollection {};
  //G4GPRObserverCollectionT<typename Type::Arg> {};

//jane fixme - do something nicer - fetch observer collection type from trigger types maybe
//template <>
//struct G4GPRObserverT<G4GPRTriggerTypes::Initialisation::RetrievePhysicsTable> : public G4GPRObserverCollectionT<G4GPRTriggerTypes::Initialisation::RetrievePhysicsTable::Arg, G4bool> {};

#endif
