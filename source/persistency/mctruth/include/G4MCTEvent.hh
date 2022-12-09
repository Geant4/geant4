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
// G4MCTEvent

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4MCTEVENT_HH
#define G4MCTEVENT_HH 1

#include <iostream>
#include <map>

#include "G4Types.hh"
#include "G4MCTGenParticle.hh"

class G4MCTGenEvent;
class G4MCTSimEvent;
class G4MCTSimParticle;

using MCTGen2SimParticleMap = std::map<G4MCTGenParticle, G4MCTSimParticle*>;
using MCTSim2GenParticleMap = std::map<G4MCTSimParticle*, G4MCTGenParticle>;

class G4MCTEvent
{
  public:

    G4MCTEvent();
    virtual ~G4MCTEvent();

    inline G4MCTEvent(const G4MCTEvent& right);
    inline G4MCTEvent& operator=(const G4MCTEvent& right);
      // copy constructor and assignment operator

    inline void SetEventNumber(G4int n);
    inline G4int GetEventNumber() const;
      // set/get functions

    inline G4MCTGenEvent* GetGenEvent() const;
    inline G4MCTSimEvent* GetSimEvent() const;

    inline G4int GetNofPrimaries() const;
    G4MCTSimParticle* GetSimParticle(const G4MCTGenParticle& genpart) const;
    G4MCTGenParticle GetGenParticle(const G4MCTSimParticle* simpart) const;
    G4int AddPrimaryPair(const G4MCTGenParticle& genp,
                        const G4MCTSimParticle* simp);
    void ClearEvent();
    void Print(std::ostream& ostr = std::cout) const;

    // iterators

    using genprimary_const_iterator = MCTGen2SimParticleMap::const_iterator;
    inline genprimary_const_iterator genprimaries_begin() const;
    inline genprimary_const_iterator genprimaries_end() const;

    using simprimary_const_iterator = MCTSim2GenParticleMap::const_iterator;
    inline simprimary_const_iterator simprimaries_begin() const;
    inline simprimary_const_iterator simprimaries_end() const;

  protected:

    G4int eventNumber = 0;
    G4MCTGenEvent* genEvent = nullptr;
    G4MCTSimEvent* simEvent = nullptr;

    // primary table (bidirectional)
    MCTGen2SimParticleMap gen2simParticleMap;
    MCTSim2GenParticleMap sim2genParticleMap;
};

// ====================================================================
// inline methods
// ====================================================================

inline G4MCTEvent::G4MCTEvent(const G4MCTEvent& right) { *this = right; }

inline G4MCTEvent& G4MCTEvent::operator=(const G4MCTEvent& right)
{
  eventNumber = right.eventNumber;

  simEvent = right.simEvent;  // shallow copy...
  genEvent = right.genEvent;

  gen2simParticleMap = right.gen2simParticleMap;
  sim2genParticleMap = right.sim2genParticleMap;

  return *this;
}

inline void G4MCTEvent::SetEventNumber(G4int n)
{
  eventNumber = n;
}

inline G4int G4MCTEvent::GetEventNumber() const
{
  return eventNumber;
}

inline G4int G4MCTEvent::GetNofPrimaries() const
{
  return (G4int)gen2simParticleMap.size();
}

inline G4MCTSimEvent* G4MCTEvent::GetSimEvent() const
{ 
  return simEvent;
}

inline G4MCTGenEvent* G4MCTEvent::GetGenEvent() const
{
  return genEvent;
}

// iterators
inline
G4MCTEvent::genprimary_const_iterator G4MCTEvent::genprimaries_begin() const
{
  return gen2simParticleMap.cbegin();
}

inline
G4MCTEvent::genprimary_const_iterator G4MCTEvent::genprimaries_end() const
{
  return gen2simParticleMap.cend();
}

inline
G4MCTEvent::simprimary_const_iterator G4MCTEvent::simprimaries_begin() const
{
  return sim2genParticleMap.cbegin();
}

inline
G4MCTEvent::simprimary_const_iterator G4MCTEvent::simprimaries_end() const
{
  return sim2genParticleMap.cend();
}

#endif
