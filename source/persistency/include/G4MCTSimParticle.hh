// $Id: G4MCTSimParticle.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// ====================================================================
//
//   G4MCTSimParticle.hh
//
// ====================================================================
#ifndef MCT_SIM_PARTICLE_H
#define MCT_SIM_PARTICLE_H

#include <vector>
#include <string>
#include <iostream>
#include "CLHEP/Vector/LorentzVector.h"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4MCTSimVertex;
class G4MCTSimParticle;

typedef std::vector<G4MCTSimParticle*> SimParticleList;

class G4MCTSimParticle {
protected:
  G4MCTSimParticle* parentParticle;
  std::vector<G4MCTSimParticle*> associatedParticleList;

  std::string name;
  int pdgID;
  int trackID;
  int parentTrackID;
  bool primaryFlag;
  HepLorentzVector momentumAtVertex;
  G4MCTSimVertex* vertex;
  bool storeFlag;
  
public:
  G4MCTSimParticle();
  G4MCTSimParticle(std::string aname, int apcode, int atid, int ptid,
		 const HepLorentzVector& p);
  G4MCTSimParticle(std::string aname, int apcode, int atid, int ptid,
		 const HepLorentzVector& p, const G4MCTSimVertex* v);
  virtual ~G4MCTSimParticle();
 
  // copy constructor and assignment operator
  G4MCTSimParticle(const G4MCTSimParticle& right);
  const G4MCTSimParticle& operator=(const G4MCTSimParticle& right);

  // set/get functions
  void SetParentParticle(const G4MCTSimParticle* p);
  G4MCTSimParticle* GetParentParticle() const;

  void SetParticleName(std::string aname);
  const std::string& GetParticleName() const;

  void SetPdgID(int id);
  int GetPdgID() const;

  void SetTrackID(int id);
  int GetTrackID() const;

  void SetParentTrackID(int id);
  int GetParentTrackID() const;

  void SetPrimaryFlag(bool q);
  bool GetPrimaryFlag() const;

  void SetMomentumAtVertex(const HepLorentzVector& p);
  const HepLorentzVector& GetMomentumAtVertex() const;

  void SetVertex(const G4MCTSimVertex* v);
  G4MCTSimVertex* GetVertex() const;

  void SetStoreFlag(bool q);
  bool GetStoreFlag() const;

  // methods...
  int AssociateParticle(G4MCTSimParticle* p);
  int GetNofAssociatedParticles() const;
  G4MCTSimParticle* GetAssociatedParticle(int i) const;
  int GetTreeLevel() const;  
  void SetStoreFlagToParentTree(bool q=true);

  void Print(std::ostream& ostr= std::cout, bool qrevorder=false) const;
  void PrintSingle(std::ostream& ostr= std::cout) const;
};

// ====================================================================
// inline functions
// ====================================================================
inline G4MCTSimParticle::G4MCTSimParticle(const G4MCTSimParticle& right)
{
  *this= right;
}
 
inline const G4MCTSimParticle& 
  G4MCTSimParticle::operator=(const G4MCTSimParticle& right)
{
  parentParticle= right.parentParticle;
  associatedParticleList= right.associatedParticleList;  // shallow copy

  name= right.name;
  pdgID= right.pdgID;
  trackID= right.trackID;
  parentTrackID= right.parentTrackID;
  primaryFlag= right.primaryFlag;
  momentumAtVertex= right.momentumAtVertex;
  vertex= right.vertex;

  return *this;
}

inline void G4MCTSimParticle::SetParentParticle(const G4MCTSimParticle* p)
{ parentParticle= const_cast<G4MCTSimParticle*>(p); }

inline G4MCTSimParticle* G4MCTSimParticle::GetParentParticle() const
{ return parentParticle; }

inline void G4MCTSimParticle::SetParticleName(std::string aname) 
{ name= aname; }

inline const std::string& G4MCTSimParticle::GetParticleName() const 
{ return name; }

inline  void G4MCTSimParticle::SetPdgID(int id) { pdgID= id; }

inline  int G4MCTSimParticle::GetPdgID() const { return pdgID; }

inline void G4MCTSimParticle::SetTrackID(int id) { trackID= id; }

inline  int G4MCTSimParticle::GetTrackID() const { return trackID; }

inline void G4MCTSimParticle::SetPrimaryFlag(bool q) { primaryFlag= q; }

inline bool G4MCTSimParticle::GetPrimaryFlag() const { return primaryFlag; }

inline void G4MCTSimParticle::SetParentTrackID(int id) 
{ parentTrackID= id; }

inline  int G4MCTSimParticle::GetParentTrackID() const 
{ return parentTrackID; }

inline  void G4MCTSimParticle::SetMomentumAtVertex(const HepLorentzVector& p)
{ momentumAtVertex= p; }

inline  const HepLorentzVector& G4MCTSimParticle::GetMomentumAtVertex() const
{ return momentumAtVertex; }

inline  void G4MCTSimParticle::SetVertex(const G4MCTSimVertex* v)
{ vertex= const_cast<G4MCTSimVertex*>(v); }

inline G4MCTSimVertex* G4MCTSimParticle::GetVertex() const
{ return vertex; }

inline void G4MCTSimParticle::SetStoreFlag(bool q) { storeFlag= q; }

inline bool G4MCTSimParticle::GetStoreFlag() const { return storeFlag; }

#endif
