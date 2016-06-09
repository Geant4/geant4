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
//   G4MCTSimParticle.hh
//
// ====================================================================
#ifndef MCT_SIM_PARTICLE_H
#define MCT_SIM_PARTICLE_H

#include "G4Types.hh"
#include <vector>
#include <string>
#include <iostream>
#include "G4LorentzVector.hh"

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
  G4bool primaryFlag;
  G4LorentzVector momentumAtVertex;
  G4MCTSimVertex* vertex;
  G4bool storeFlag;
  
public:
  G4MCTSimParticle();
  G4MCTSimParticle(std::string aname, int apcode, int atid, int ptid,
		 const G4LorentzVector& p);
  G4MCTSimParticle(std::string aname, int apcode, int atid, int ptid,
		 const G4LorentzVector& p, const G4MCTSimVertex* v);
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

  void SetPrimaryFlag(G4bool q);
  G4bool GetPrimaryFlag() const;

  void SetMomentumAtVertex(const G4LorentzVector& p);
  const G4LorentzVector& GetMomentumAtVertex() const;

  void SetVertex(const G4MCTSimVertex* v);
  G4MCTSimVertex* GetVertex() const;

  void SetStoreFlag(G4bool q);
  G4bool GetStoreFlag() const;

  // methods...
  int AssociateParticle(G4MCTSimParticle* p);
  int GetNofAssociatedParticles() const;
  G4MCTSimParticle* GetAssociatedParticle(int i) const;
  int GetTreeLevel() const;  
  void SetStoreFlagToParentTree(G4bool q=true);

  void Print(std::ostream& ostr= std::cout, G4bool qrevorder=false) const;
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

inline void G4MCTSimParticle::SetPrimaryFlag(G4bool q) { primaryFlag= q; }

inline G4bool G4MCTSimParticle::GetPrimaryFlag() const { return primaryFlag; }

inline void G4MCTSimParticle::SetParentTrackID(int id) 
{ parentTrackID= id; }

inline  int G4MCTSimParticle::GetParentTrackID() const 
{ return parentTrackID; }

inline  void G4MCTSimParticle::SetMomentumAtVertex(const G4LorentzVector& p)
{ momentumAtVertex= p; }

inline  const G4LorentzVector& G4MCTSimParticle::GetMomentumAtVertex() const
{ return momentumAtVertex; }

inline  void G4MCTSimParticle::SetVertex(const G4MCTSimVertex* v)
{ vertex= const_cast<G4MCTSimVertex*>(v); }

inline G4MCTSimVertex* G4MCTSimParticle::GetVertex() const
{ return vertex; }

inline void G4MCTSimParticle::SetStoreFlag(G4bool q) { storeFlag= q; }

inline G4bool G4MCTSimParticle::GetStoreFlag() const { return storeFlag; }

#endif
