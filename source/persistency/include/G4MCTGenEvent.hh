// $Id: G4MCTGenEvent.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// ====================================================================
//
//   G4MCTGenEvent.hh
//
// ====================================================================
#ifndef MCT_GEN_EVENT_H
#define MCT_GEN_EVENT_H

#include "G4Types.hh"
#include "g4std/iostream"
#include "g4std/vector"
#include "CLHEP/HepMC/GenEvent.h"
 
// ====================================================================
//
// class definition
//
// ====================================================================

class G4MCTGenEvent {
protected:
  G4std::vector<HepMC::GenEvent*> eventList;

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

  void Print(G4std::ostream& ostr= G4std::cout) const;
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
