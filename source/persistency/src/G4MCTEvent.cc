//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//   G4MCTEvent.cc
//
// ====================================================================

#ifndef WIN32

#include "G4ios.hh"
#include "G4MCTEvent.hh"
#include "G4MCTGenEvent.hh"
#include "G4MCTSimEvent.hh"
#include "G4MCTSimParticle.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////
G4MCTEvent::G4MCTEvent()
  : eventNumber(0)
////////////////////
{
  genEvent= new G4MCTGenEvent();  // will be reused.
  simEvent= new G4MCTSimEvent();
}

/////////////////////
G4MCTEvent::~G4MCTEvent()
/////////////////////
{
  delete genEvent;
  delete simEvent;
}

/////////////////////////////////////////////////////
G4MCTSimParticle* G4MCTEvent::GetSimParticle
                (const G4MCTGenParticle& genpart) const
/////////////////////////////////////////////////////
{
  MCTGen2SimParticleMap::const_iterator pos= gen2simParticleMap.find(genpart);
  if(pos != gen2simParticleMap.end()) {
    return pos-> second;
  } else {
    return 0;
  }
}

////////////////////////////////////////////////////
G4MCTGenParticle G4MCTEvent::GetGenParticle
               (const G4MCTSimParticle* simpart) const
////////////////////////////////////////////////////
{
  MCTSim2GenParticleMap::const_iterator 
    pos= sim2genParticleMap.find(const_cast<G4MCTSimParticle*>(simpart));
  if(pos != sim2genParticleMap.end()) {
    return pos-> second;
  } else {
    return G4MCTGenParticle(0,0);
  }
}

////////////////////////////////////////////////////////
int G4MCTEvent::AddPrimaryPair(const G4MCTGenParticle& genp, 
			     const G4MCTSimParticle* simp)
////////////////////////////////////////////////////////
{
  gen2simParticleMap.insert(G4std::make_pair(const_cast<G4MCTGenParticle&>(genp), 
					   const_cast<G4MCTSimParticle*>(simp)));
  sim2genParticleMap.insert(G4std::make_pair(const_cast<G4MCTSimParticle*>(simp), 
					   const_cast<G4MCTGenParticle&>(genp)));

  return gen2simParticleMap.size();
}

///////////////////////////
void G4MCTEvent::ClearEvent()
///////////////////////////
{
  gen2simParticleMap.clear();
  sim2genParticleMap.clear();

  genEvent-> ClearEvent();
  simEvent-> ClearEvent();
}


//////////////////////////////////////////////
void G4MCTEvent::Print(G4std::ostream& ostr) const
//////////////////////////////////////////////
{
  ostr << "Event#:" << eventNumber << G4endl;
  genEvent-> Print(ostr);
  simEvent-> Print(ostr);
}

#endif
