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
// G4MCTSimEvent

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4MCTSIMEVENT_HH
#define G4MCTSIMEVENT_HH 1

#include <iostream>
#include <vector>
#include <map>

#include "G4Types.hh"

class G4MCTSimParticle;
class G4MCTSimVertex;

using G4MCTSimParticleContainer = std::map<int, G4MCTSimParticle*>;
using G4MCTSimVertexContainer = std::vector<G4MCTSimVertex*>;

class G4MCTSimEvent
{
  public:

    G4MCTSimEvent();
    ~G4MCTSimEvent();

    inline G4MCTSimEvent(const G4MCTSimEvent& right);
    inline G4MCTSimEvent& operator=(const G4MCTSimEvent& right);
      // copy constructor and assignment operator

    G4bool AddParticle(const G4MCTSimParticle* aparticle);
    inline G4int GetNofParticles() const;
    inline G4int GetNofVertices() const;
    G4int GetNofStoredParticles() const;
    G4int GetNofStoredVertices() const;
    G4MCTSimParticle* FindParticle(G4int tid) const;
    G4MCTSimVertex* GetVertex(G4int vid) const;

    void BuildVertexContainer();
    void ClearEvent();
    void Print(std::ostream& ostr = std::cout) const;

    // iterators
    using particle_iterator = G4MCTSimParticleContainer::iterator;
    using particle_const_iterator = G4MCTSimParticleContainer::const_iterator;
    inline particle_iterator particles_begin();
    inline particle_iterator particles_end();
    inline particle_const_iterator particles_begin() const;
    inline particle_const_iterator particles_end() const;

    using vertex_iterator = G4MCTSimVertexContainer::iterator;
    using vertex_const_iterator = G4MCTSimVertexContainer::const_iterator;
    inline vertex_iterator vertices_begin();
    inline vertex_iterator vertices_end();
    inline vertex_const_iterator vertices_begin() const;
    inline vertex_const_iterator vertices_end() const;

  protected:

    G4MCTSimParticleContainer particleMap;
    G4MCTSimVertexContainer vertexVec;
};

// ====================================================================
// inline methods
// ====================================================================

inline G4MCTSimEvent::G4MCTSimEvent(const G4MCTSimEvent& right)
{
  *this = right;
}

inline G4MCTSimEvent& G4MCTSimEvent::operator=(const G4MCTSimEvent& right)
{
  particleMap = right.particleMap;  // shallow copy

  return *this;
}

inline G4int G4MCTSimEvent::GetNofParticles() const
{
  return (G4int)particleMap.size();
}

inline G4int G4MCTSimEvent::GetNofVertices() const
{
  return (G4int)vertexVec.size();
}

// iterators
inline G4MCTSimEvent::particle_iterator G4MCTSimEvent::particles_begin()
{
  return particleMap.begin();
}

inline G4MCTSimEvent::particle_iterator G4MCTSimEvent::particles_end()
{
  return particleMap.end();
}

inline G4MCTSimEvent::particle_const_iterator G4MCTSimEvent::particles_begin()
  const
{
  return particleMap.cbegin();
}

inline G4MCTSimEvent::particle_const_iterator G4MCTSimEvent::particles_end()
  const
{
  return particleMap.cend();
}

inline G4MCTSimEvent::vertex_iterator G4MCTSimEvent::vertices_begin()
{
  return vertexVec.begin();
}

inline G4MCTSimEvent::vertex_iterator G4MCTSimEvent::vertices_end()
{
  return vertexVec.end();
}

inline G4MCTSimEvent::vertex_const_iterator G4MCTSimEvent::vertices_begin()
  const
{
  return vertexVec.cbegin();
}

inline G4MCTSimEvent::vertex_const_iterator G4MCTSimEvent::vertices_end() const
{
  return vertexVec.cend();
}

#endif
