// $Id: G4MCTEvent.cc,v 1.1 2002-11-24 13:45:24 morita Exp $
// ====================================================================
//
//   G4MCTEvent.cc
//
// ====================================================================
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
  gen2simParticleMap.insert(std::make_pair(const_cast<G4MCTGenParticle&>(genp), 
					   const_cast<G4MCTSimParticle*>(simp)));
  sim2genParticleMap.insert(std::make_pair(const_cast<G4MCTSimParticle*>(simp), 
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
void G4MCTEvent::Print(std::ostream& ostr) const
//////////////////////////////////////////////
{
  ostr << "Event#:" << eventNumber << std::endl;
  genEvent-> Print(ostr);
  simEvent-> Print(ostr);
}
