// $Id: G4MCTSimVertex.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// ====================================================================
//
//   G4MCTSimVertex.hh
//
// ====================================================================
#ifndef MCT_SIM_VERTEX_H
#define MCT_SIM_VERTEX_H

#include <iostream>
#include <vector>
#include <string>
#include "CLHEP/Vector/ThreeVector.h"
#include "G4MCTSimParticle.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class G4MCTSimVertex {
private:
  int  inParticleTrackID;
  std::vector<int> outParticleTrackIDList;

  int id; // assigned independently from G4
  Hep3Vector position;
  double time;
  std::string volumeName;
  int volumeNumber;
  std::string creatorProcessName;
  bool storeFlag;
  
public:
  G4MCTSimVertex();
  G4MCTSimVertex(const Hep3Vector& x, double t);
  G4MCTSimVertex(const Hep3Vector& x, double t, 
	       std::string vname, int ncopy, std::string pname);
  ~G4MCTSimVertex();
 
  // copy constructor and assignment operator
  G4MCTSimVertex(const G4MCTSimVertex& right);
  const G4MCTSimVertex& operator=(const G4MCTSimVertex& right);

  // set/get functions
  void SetID(int i);
  int GetID() const;

  void SetPosition(const Hep3Vector& x);
  const Hep3Vector& GetPosition() const;

  void SetTime(double t);
  double GetTime() const;

  void SetVolumeName(std::string vname);
  const std::string& GetVolumeName() const;

  void SetVolumeNumber(int n);
  int GetVolumeNumber() const;

  void SetCreatorProcessName(std::string pname);
  const std::string& GetCreatorProcessName() const;

  void SetStoreFlag(bool q);
  bool GetStoreFlag() const;

  // methods...
  void SetInParticle(const G4MCTSimParticle* in);
  void SetInParticle(int in);
  int GetInParticleTrackID() const;

  int GetNofOutParticles() const;
  int AddOutParticle(const G4MCTSimParticle* out);
  int AddOutParticle(int out);
  int GetOutParticleTrackID(int i) const;

  void Print(std::ostream& ostr= std::cout) const;
};

// ====================================================================
// inline functions
// ====================================================================
inline G4MCTSimVertex::G4MCTSimVertex(const G4MCTSimVertex& right)
{
  *this= right;
}
 
inline const G4MCTSimVertex& 
  G4MCTSimVertex::operator=(const G4MCTSimVertex& right)
{
  inParticleTrackID= right.inParticleTrackID;
  outParticleTrackIDList= right.outParticleTrackIDList;

  id= right.id;
  position= right.position;
  time= right.time;
  volumeName= right.volumeName;
  volumeNumber= right.volumeNumber;
  creatorProcessName= right.creatorProcessName;

  return *this;
}

inline void G4MCTSimVertex::SetID(int i) { id= i; }
inline int  G4MCTSimVertex::GetID() const { return id; }

inline  void G4MCTSimVertex::SetPosition(const Hep3Vector& x)
{ position= x; }

inline  const Hep3Vector& G4MCTSimVertex::GetPosition() const
{ return position; }

inline  void G4MCTSimVertex::SetTime(double t)
{ time= t; }

inline double G4MCTSimVertex::GetTime() const
{ return time; }

inline  void G4MCTSimVertex::SetVolumeName(std::string vname)
{ volumeName= vname; }

inline  const std::string& G4MCTSimVertex::GetVolumeName() const
{ return volumeName; }

inline void G4MCTSimVertex::SetVolumeNumber(int n)
{ volumeNumber= n; }

inline int G4MCTSimVertex::GetVolumeNumber() const
{ return volumeNumber; }

inline  void G4MCTSimVertex::SetCreatorProcessName(std::string pname)
{ creatorProcessName= pname; }

inline  const std::string& G4MCTSimVertex::GetCreatorProcessName() const
{ return creatorProcessName; }

inline void G4MCTSimVertex::SetStoreFlag(bool q) { storeFlag= q; }

inline bool G4MCTSimVertex::GetStoreFlag() const { return storeFlag; }

inline void G4MCTSimVertex::SetInParticle(const G4MCTSimParticle* in)
{ inParticleTrackID= in-> GetTrackID(); }

inline void G4MCTSimVertex::SetInParticle(int in)
{ inParticleTrackID= in; }

inline int G4MCTSimVertex::GetInParticleTrackID() const
{ return inParticleTrackID; }

inline int G4MCTSimVertex::GetNofOutParticles() const
{ return outParticleTrackIDList.size(); }

inline int G4MCTSimVertex::AddOutParticle(const G4MCTSimParticle* out)
{ 
  outParticleTrackIDList.push_back(out->GetTrackID());
  return outParticleTrackIDList.size();
}

inline int G4MCTSimVertex::AddOutParticle(int out)
{
  outParticleTrackIDList.push_back(out);
  return outParticleTrackIDList.size();
}

inline int G4MCTSimVertex::GetOutParticleTrackID(int i) const
{
  int size= outParticleTrackIDList.size();
  if(i>=0 && i< size) return outParticleTrackIDList[i];
  else return 0;    
}

#endif
