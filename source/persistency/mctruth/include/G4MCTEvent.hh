//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//   G4MCTEvent.hh
//
// ====================================================================
#ifndef MCT_EVENT_H
#define MCT_EVENT_H

#include "G4Types.hh"
#include <iostream>
#include <map>
#include "G4MCTGenParticle.hh"
 
// ====================================================================
//
// class definition
//
// ====================================================================
class G4MCTGenEvent;
class G4MCTSimEvent;
class G4MCTSimParticle;

typedef std::map<G4MCTGenParticle, G4MCTSimParticle*> MCTGen2SimParticleMap;
typedef std::map<G4MCTSimParticle*, G4MCTGenParticle> MCTSim2GenParticleMap;

class G4MCTEvent {
protected:
  int eventNumber;
  G4MCTGenEvent* genEvent;
  G4MCTSimEvent* simEvent;

  // primary table (bidirectional)
  MCTGen2SimParticleMap gen2simParticleMap;
  MCTSim2GenParticleMap sim2genParticleMap;

public:
  G4MCTEvent();
  virtual ~G4MCTEvent();
 
  // copy constructor and assignment operator
  G4MCTEvent(const G4MCTEvent& right);
  const G4MCTEvent& operator=(const G4MCTEvent& right);

  // set/get functions
  void SetEventNumber(int n);
  int GetEventNumber() const;

  G4MCTGenEvent* GetGenEvent() const;
  G4MCTSimEvent* GetSimEvent() const;

  // methods...
  int GetNofPrimaries() const;
  G4MCTSimParticle* GetSimParticle(const G4MCTGenParticle& genpart) const;
  G4MCTGenParticle GetGenParticle(const G4MCTSimParticle* simpart) const;
  int AddPrimaryPair(const G4MCTGenParticle& genp, 
		     const G4MCTSimParticle* simp); 
  void ClearEvent();
  void Print(std::ostream& ostr= std::cout) const;

  // iterators
  typedef MCTGen2SimParticleMap::const_iterator genprimary_const_iterator;
  genprimary_const_iterator genprimaries_begin() const;
  genprimary_const_iterator genprimaries_end() const;
  
  typedef MCTSim2GenParticleMap::const_iterator simprimary_const_iterator;
  simprimary_const_iterator simprimaries_begin() const;
  simprimary_const_iterator simprimaries_end() const;
};

// ====================================================================
// inline functions
// ====================================================================

inline G4MCTEvent::G4MCTEvent(const G4MCTEvent& right)
{
  *this= right;
}
 
inline const G4MCTEvent& G4MCTEvent::operator=(const G4MCTEvent& right)
{
  eventNumber= right.eventNumber;

  simEvent= right.simEvent; // shallow copy...
  genEvent= right.genEvent;

  gen2simParticleMap= right.gen2simParticleMap;
  sim2genParticleMap= right.sim2genParticleMap;

  return *this;
}

inline void G4MCTEvent::SetEventNumber(int n) { eventNumber= n; }
inline int G4MCTEvent::GetEventNumber() const { return eventNumber; }

inline int G4MCTEvent::GetNofPrimaries() const 
           { return gen2simParticleMap.size(); }
inline G4MCTSimEvent* G4MCTEvent::GetSimEvent() const { return simEvent; }
inline G4MCTGenEvent* G4MCTEvent::GetGenEvent() const { return genEvent; }

// iterators
inline G4MCTEvent::genprimary_const_iterator G4MCTEvent::genprimaries_begin() const
{ return gen2simParticleMap.begin(); }

inline G4MCTEvent::genprimary_const_iterator G4MCTEvent::genprimaries_end() const
{ return gen2simParticleMap.end(); }

inline G4MCTEvent::simprimary_const_iterator G4MCTEvent::simprimaries_begin() const
{ return sim2genParticleMap.begin(); }

inline G4MCTEvent::simprimary_const_iterator G4MCTEvent::simprimaries_end() const
{ return sim2genParticleMap.end(); }

#endif
