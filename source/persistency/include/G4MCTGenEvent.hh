// $Id: G4MCTGenEvent.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// ====================================================================
//
//   G4MCTGenEvent.hh
//
// ====================================================================
#ifndef MCT_GEN_EVENT_H
#define MCT_GEN_EVENT_H

#include <iostream>
#include <vector>
#include "CLHEP/HepMC/GenEvent.h"
 
// ====================================================================
//
// class definition
//
// ====================================================================

class G4MCTGenEvent {
protected:
  std::vector<HepMC::GenEvent*> eventList;

public:
  G4MCTGenEvent();
  virtual ~G4MCTGenEvent();
 
  // copy constructor and assignment operator
  G4MCTGenEvent(const G4MCTGenEvent& right);
  const G4MCTGenEvent& operator=(const G4MCTGenEvent& right);
  
  // methods...  
  int AddGenEvent(const HepMC::GenEvent* genevent);
  int GetNofEvents() const;
  const HepMC::GenEvent* GetGenEvent(int i);

  void ClearEvent();

  void Print(std::ostream& ostr= std::cout) const;
};

// ====================================================================
// inline functions
// ====================================================================

inline G4MCTGenEvent::G4MCTGenEvent(const G4MCTGenEvent& right)
{
  *this= right;
}
 
inline const G4MCTGenEvent& G4MCTGenEvent::operator=(const G4MCTGenEvent& right)
{
  eventList= right.eventList;  // shallow copy

  return *this;
}

#endif
