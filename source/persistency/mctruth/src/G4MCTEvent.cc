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
// G4MCTEvent implementation
//
// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------

#include "G4ios.hh"
#include "G4MCTEvent.hh"
#include "G4MCTGenEvent.hh"
#include "G4MCTSimEvent.hh"
#include "G4MCTSimParticle.hh"

// --------------------------------------------------------------------
G4MCTEvent::G4MCTEvent()
{
  genEvent = new G4MCTGenEvent();  // will be reused.
  simEvent = new G4MCTSimEvent();
}

// --------------------------------------------------------------------
G4MCTEvent::~G4MCTEvent()
{
  delete genEvent;
  delete simEvent;
}

// --------------------------------------------------------------------
G4MCTSimParticle* G4MCTEvent::GetSimParticle(
  const G4MCTGenParticle& genpart) const
{
  auto pos = gen2simParticleMap.find(genpart);
  if(pos != gen2simParticleMap.cend())
  {
    return pos->second;
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
G4MCTGenParticle G4MCTEvent::GetGenParticle(
  const G4MCTSimParticle* simpart) const
{
  auto pos = sim2genParticleMap.find(const_cast<G4MCTSimParticle*>(simpart));
  if(pos != sim2genParticleMap.cend())
  {
    return pos->second;
  }
  else
  {
    return G4MCTGenParticle((void*) 0, (void*) 0);
  }
}

// --------------------------------------------------------------------
G4int G4MCTEvent::AddPrimaryPair(const G4MCTGenParticle& genp,
                                 const G4MCTSimParticle* simp)
{
  gen2simParticleMap.insert(std::make_pair(
    const_cast<G4MCTGenParticle&>(genp), const_cast<G4MCTSimParticle*>(simp)));
  sim2genParticleMap.insert(std::make_pair(
    const_cast<G4MCTSimParticle*>(simp), const_cast<G4MCTGenParticle&>(genp)));

  return (G4int)gen2simParticleMap.size();
}

// --------------------------------------------------------------------
void G4MCTEvent::ClearEvent()
{
  gen2simParticleMap.clear();
  sim2genParticleMap.clear();

  genEvent->ClearEvent();
  simEvent->ClearEvent();
}

// --------------------------------------------------------------------
void G4MCTEvent::Print(std::ostream& ostr) const
{
  ostr << "Event#:" << eventNumber << G4endl;
  simEvent->Print(ostr);
}
