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
//   G4MCTSimEvent.hh
//
// ====================================================================
#ifndef MCT_SIM_EVENT_H
#define MCT_SIM_EVENT_H

#include "G4Types.hh"
#include "g4std/iostream"
#include "g4std/vector"
#include "g4std/map"
 
// ====================================================================
//
// class definition
//
// ====================================================================
class G4MCTSimParticle;
class G4MCTSimVertex;

typedef G4std::map<int, G4MCTSimParticle*> G4MCTSimParticleContainer;
typedef G4std::vector<G4MCTSimVertex*> G4MCTSimVertexContainer;

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
  void Print(G4std::ostream& ostr= G4std::cout) const;  

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
