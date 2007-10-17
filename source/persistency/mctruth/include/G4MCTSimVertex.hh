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
//   G4MCTSimVertex.hh
//
// ====================================================================
#ifndef MCT_SIM_VERTEX_H
#define MCT_SIM_VERTEX_H

#include "G4Types.hh"
#include <iostream>
#include <vector>
#include <string>
#include "G4ThreeVector.hh"
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
  G4ThreeVector position;
  double time;
  std::string volumeName;
  int volumeNumber;
  std::string creatorProcessName;
  G4bool storeFlag;
  
public:
  G4MCTSimVertex();
  G4MCTSimVertex(const G4ThreeVector& x, double t);
  G4MCTSimVertex(const G4ThreeVector& x, double t, 
	       std::string vname, int ncopy, std::string pname);
  ~G4MCTSimVertex();
 
  // copy constructor and assignment operator
  G4MCTSimVertex(const G4MCTSimVertex& right);
  const G4MCTSimVertex& operator=(const G4MCTSimVertex& right);

  // set/get functions
  void SetID(int i);
  int GetID() const;

  void SetPosition(const G4ThreeVector& x);
  const G4ThreeVector& GetPosition() const;

  void SetTime(double t);
  double GetTime() const;

  void SetVolumeName(std::string vname);
  const std::string& GetVolumeName() const;

  void SetVolumeNumber(int n);
  int GetVolumeNumber() const;

  void SetCreatorProcessName(std::string pname);
  const std::string& GetCreatorProcessName() const;

  void SetStoreFlag(G4bool q);
  G4bool GetStoreFlag() const;

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

inline  void G4MCTSimVertex::SetPosition(const G4ThreeVector& x)
{ position= x; }

inline  const G4ThreeVector& G4MCTSimVertex::GetPosition() const
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

inline void G4MCTSimVertex::SetStoreFlag(G4bool q) { storeFlag= q; }

inline G4bool G4MCTSimVertex::GetStoreFlag() const { return storeFlag; }

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
