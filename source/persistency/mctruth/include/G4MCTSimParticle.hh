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
// G4MCTSimParticle

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4MCTSIMPARTICLE_HH
#define G4MCTSIMPARTICLE_HH 1

#include <vector>
#include <string>
#include <iostream>

#include "G4String.hh"
#include "G4Types.hh"
#include "G4LorentzVector.hh"

class G4MCTSimVertex;
class G4MCTSimParticle;

using SimParticleList = std::vector<G4MCTSimParticle*>;

class G4MCTSimParticle
{
  public:

    G4MCTSimParticle();
    G4MCTSimParticle(const G4String& aname,
                     G4int apcode, G4int atid, G4int ptid,
                     const G4LorentzVector& p);
    G4MCTSimParticle(const G4String& aname,
                     G4int apcode, G4int atid, G4int ptid,
                     const G4LorentzVector& p, const G4MCTSimVertex* v);
    virtual ~G4MCTSimParticle();

    inline G4MCTSimParticle(const G4MCTSimParticle& right);
    inline G4MCTSimParticle& operator=(const G4MCTSimParticle& right);
      // copy constructor and assignment operator

    inline void SetParentParticle(const G4MCTSimParticle* p);
    inline G4MCTSimParticle* GetParentParticle() const;

    inline void SetParticleName(std::string aname);
    inline const G4String& GetParticleName() const;

    inline void SetPdgID(G4int id);
    inline G4int GetPdgID() const;

    inline void SetTrackID(G4int id);
    inline G4int GetTrackID() const;

    inline void SetParentTrackID(G4int id);
    inline G4int GetParentTrackID() const;

    inline void SetPrimaryFlag(G4bool q);
    inline G4bool GetPrimaryFlag() const;

    inline void SetMomentumAtVertex(const G4LorentzVector& p);
    inline const G4LorentzVector& GetMomentumAtVertex() const;

    inline void SetVertex(const G4MCTSimVertex* v);
    inline G4MCTSimVertex* GetVertex() const;

    inline void SetStoreFlag(G4bool q);
    inline G4bool GetStoreFlag() const;

    G4int AssociateParticle(G4MCTSimParticle* p);
    G4int GetNofAssociatedParticles() const;
    G4MCTSimParticle* GetAssociatedParticle(G4int i) const;
    G4int GetTreeLevel() const;
    void SetStoreFlagToParentTree(G4bool q = true);

    void Print(std::ostream& ostr = std::cout, G4bool qrevorder = false) const;
    void PrintSingle(std::ostream& ostr = std::cout) const;

  protected:

    G4MCTSimParticle* parentParticle = nullptr;
    std::vector<G4MCTSimParticle*> associatedParticleList;

    G4String name;
    G4LorentzVector momentumAtVertex;
    G4MCTSimVertex* vertex = nullptr;
    G4int pdgID = 0;
    G4int trackID = 0;
    G4int parentTrackID = 0;
    G4bool primaryFlag = false;
    G4bool storeFlag = false;
};

// ====================================================================
// inline methods
// ====================================================================

inline G4MCTSimParticle::G4MCTSimParticle(const G4MCTSimParticle& right)
{
  *this = right;
}

inline G4MCTSimParticle& G4MCTSimParticle::operator=(
  const G4MCTSimParticle& right)
{
  parentParticle         = right.parentParticle;
  associatedParticleList = right.associatedParticleList;  // shallow copy

  name             = right.name;
  pdgID            = right.pdgID;
  trackID          = right.trackID;
  parentTrackID    = right.parentTrackID;
  primaryFlag      = right.primaryFlag;
  momentumAtVertex = right.momentumAtVertex;
  vertex           = right.vertex;

  return *this;
}

inline void G4MCTSimParticle::SetParentParticle(const G4MCTSimParticle* p)
{
  parentParticle = const_cast<G4MCTSimParticle*>(p);
}

inline G4MCTSimParticle* G4MCTSimParticle::GetParentParticle() const
{
  return parentParticle;
}

inline void G4MCTSimParticle::SetParticleName(std::string aname)
{
  name = aname;
}

inline const G4String& G4MCTSimParticle::GetParticleName() const
{
  return name;
}

inline void G4MCTSimParticle::SetPdgID(G4int id)
{
  pdgID = id;
}

inline G4int G4MCTSimParticle::GetPdgID() const
{
  return pdgID;
}

inline void G4MCTSimParticle::SetTrackID(G4int id)
{
  trackID = id;
}

inline G4int G4MCTSimParticle::GetTrackID() const
{
  return trackID;
}

inline void G4MCTSimParticle::SetPrimaryFlag(G4bool q)
{
  primaryFlag = q;
}

inline G4bool G4MCTSimParticle::GetPrimaryFlag() const
{
  return primaryFlag;
}

inline void G4MCTSimParticle::SetParentTrackID(G4int id)
{
  parentTrackID = id;
}

inline G4int G4MCTSimParticle::GetParentTrackID() const
{
  return parentTrackID;
}

inline void G4MCTSimParticle::SetMomentumAtVertex(const G4LorentzVector& p)
{
  momentumAtVertex = p;
}

inline const G4LorentzVector& G4MCTSimParticle::GetMomentumAtVertex() const
{
  return momentumAtVertex;
}

inline void G4MCTSimParticle::SetVertex(const G4MCTSimVertex* v)
{
  vertex = const_cast<G4MCTSimVertex*>(v);
}

inline G4MCTSimVertex* G4MCTSimParticle::GetVertex() const
{
  return vertex;
}

inline void G4MCTSimParticle::SetStoreFlag(G4bool q)
{
  storeFlag = q;
}

inline G4bool G4MCTSimParticle::GetStoreFlag() const
{
  return storeFlag;
}

#endif
