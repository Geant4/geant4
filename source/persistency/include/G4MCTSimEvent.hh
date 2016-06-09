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
//   G4MCTSimEvent.hh
//
// ====================================================================
#ifndef MCT_SIM_EVENT_H
#define MCT_SIM_EVENT_H

#include "G4Types.hh"
#include <iostream>
#include <vector>
#include <map>
 
// ====================================================================
//
// class definition
//
// ====================================================================
class G4MCTSimParticle;
class G4MCTSimVertex;

typedef std::map<int, G4MCTSimParticle*> G4MCTSimParticleContainer;
typedef std::vector<G4MCTSimVertex*> G4MCTSimVertexContainer;

class G4MCTSimEvent {
protected:
G4MCTSimParticleContainer particleMap;
G4MCTSimVertexContainer vertexVec;
 
public:
  G4MCTSimEvent();
  ~G4MCTSimEvent();
 
  // copy constructor and assignment operator
  G4MCTSimEvent(const G4MCTSimEvent& right);
  const G4MCTSimEvent& operator=(const G4MCTSimEvent& right);

  // methods...
  G4bool AddParticle(const G4MCTSimParticle* aparticle);
  int GetNofParticles() const;
  int GetNofVertices() const;
  int GetNofStoredParticles() const;
  int GetNofStoredVertices() const;
  G4MCTSimParticle* FindParticle(int tid) const;
  G4MCTSimVertex* GetVertex(int vid) const;

  void BuildVertexContainer();
  void ClearEvent();
  void Print(std::ostream& ostr= std::cout) const;  

  // iterators
  typedef G4MCTSimParticleContainer::iterator particle_iterator;
  typedef G4MCTSimParticleContainer::const_iterator particle_const_iterator;
  particle_iterator particles_begin();
  particle_iterator particles_end();
  particle_const_iterator particles_begin() const;
  particle_const_iterator particles_end() const;
  
  typedef G4MCTSimVertexContainer::iterator vertex_iterator;
  typedef G4MCTSimVertexContainer::const_iterator vertex_const_iterator;
  vertex_iterator vertices_begin();
  vertex_iterator vertices_end();
  vertex_const_iterator vertices_begin() const;
  vertex_const_iterator vertices_end() const;
};

// ====================================================================
// inline functions
// ====================================================================

inline G4MCTSimEvent::G4MCTSimEvent(const G4MCTSimEvent& right)
{
  *this= right;
}
 
inline const G4MCTSimEvent& G4MCTSimEvent::operator=(const G4MCTSimEvent& right)
{
  particleMap= right.particleMap; // shallow copy

  return *this;
}

inline int G4MCTSimEvent::GetNofParticles() const
{
  return particleMap.size();
}

inline int G4MCTSimEvent::GetNofVertices() const
{
  return vertexVec.size();
}

// iterators
inline G4MCTSimEvent::particle_iterator G4MCTSimEvent::particles_begin()
{ return particleMap.begin(); }

inline G4MCTSimEvent::particle_iterator G4MCTSimEvent::particles_end()
{ return particleMap.end(); }

inline G4MCTSimEvent::particle_const_iterator 
                    G4MCTSimEvent::particles_begin() const
{ return particleMap.begin(); }

inline G4MCTSimEvent::particle_const_iterator G4MCTSimEvent::particles_end() const
{ return particleMap.end(); }

inline G4MCTSimEvent::vertex_iterator G4MCTSimEvent::vertices_begin()
{ return vertexVec.begin(); }

inline G4MCTSimEvent::vertex_iterator G4MCTSimEvent::vertices_end()
{ return vertexVec.end(); }

inline G4MCTSimEvent::vertex_const_iterator G4MCTSimEvent::vertices_begin() const
{ return vertexVec.begin(); }

inline G4MCTSimEvent::vertex_const_iterator G4MCTSimEvent::vertices_end() const
{ return vertexVec.end(); }

#endif
