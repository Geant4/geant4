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
// G4MCTSimVertex

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4MCTSIMVERTEX_HH
#define G4MCTSIMVERTEX_HH 1

#include <iostream>
#include <vector>
#include <string>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MCTSimParticle.hh"

class G4MCTSimVertex
{
  public:

    G4MCTSimVertex();
    G4MCTSimVertex(const G4ThreeVector& x, G4double t);
    G4MCTSimVertex(const G4ThreeVector& x, G4double t,
                   const G4String& vname, G4int ncopy,
                   const G4String& pname);
    ~G4MCTSimVertex();

    inline G4MCTSimVertex(const G4MCTSimVertex& right);
    inline G4MCTSimVertex& operator=(const G4MCTSimVertex& right);
      // copy constructor and assignment operator

    inline void SetID(G4int i);
    inline G4int GetID() const;

    inline void SetPosition(const G4ThreeVector& x);
    inline const G4ThreeVector& GetPosition() const;

    inline void SetTime(G4double t);
    inline G4double GetTime() const;

    inline void SetVolumeName(const G4String& vname);
    inline const G4String& GetVolumeName() const;

    inline void SetVolumeNumber(G4int n);
    inline G4int GetVolumeNumber() const;

    inline void SetCreatorProcessName(const G4String& pname);
    inline const G4String& GetCreatorProcessName() const;

    inline void SetStoreFlag(G4bool q);
    inline G4bool GetStoreFlag() const;

    inline void SetInParticle(const G4MCTSimParticle* in);
    inline void SetInParticle(G4int in);
    inline G4int GetInParticleTrackID() const;

    inline G4int GetNofOutParticles() const;
    inline G4int AddOutParticle(const G4MCTSimParticle* out);
    inline G4int AddOutParticle(G4int out);
    inline G4int GetOutParticleTrackID(G4int i) const;

    void Print(std::ostream& ostr = std::cout) const;

  private:

    G4int inParticleTrackID = 0;
    std::vector<G4int> outParticleTrackIDList;

    G4String volumeName = "";
    G4String creatorProcessName = "none";
    G4ThreeVector position;
    G4double time = 0.0;
    G4int id = -1;  // assigned independently from G4
    G4int volumeNumber = -1;
    G4bool storeFlag = false;
};

// ====================================================================
// inline methods
// ====================================================================

inline G4MCTSimVertex::G4MCTSimVertex(const G4MCTSimVertex& right)
{
  *this = right;
}

inline G4MCTSimVertex& G4MCTSimVertex::operator=(
  const G4MCTSimVertex& right)
{
  inParticleTrackID      = right.inParticleTrackID;
  outParticleTrackIDList = right.outParticleTrackIDList;

  id                 = right.id;
  position           = right.position;
  time               = right.time;
  volumeName         = right.volumeName;
  volumeNumber       = right.volumeNumber;
  creatorProcessName = right.creatorProcessName;

  return *this;
}

inline void G4MCTSimVertex::SetID(G4int i)
{
  id = i;
}

inline G4int G4MCTSimVertex::GetID() const
{
  return id;
}

inline void G4MCTSimVertex::SetPosition(const G4ThreeVector& x)
{
  position = x;
}

inline const G4ThreeVector& G4MCTSimVertex::GetPosition() const
{
  return position;
}

inline void G4MCTSimVertex::SetTime(G4double t)
{
  time = t;
}

inline G4double G4MCTSimVertex::GetTime() const
{
  return time;
}

inline void G4MCTSimVertex::SetVolumeName(const G4String& vname)
{
  volumeName = vname;
}

inline const G4String& G4MCTSimVertex::GetVolumeName() const
{
  return volumeName;
}

inline void G4MCTSimVertex::SetVolumeNumber(G4int n)
{
  volumeNumber = n;
}

inline G4int G4MCTSimVertex::GetVolumeNumber() const
{
  return volumeNumber;
}

inline void G4MCTSimVertex::SetCreatorProcessName(const G4String& pname)
{
  creatorProcessName = pname;
}

inline const G4String& G4MCTSimVertex::GetCreatorProcessName() const
{
  return creatorProcessName;
}

inline void G4MCTSimVertex::SetStoreFlag(G4bool q)
{
  storeFlag = q;
}

inline G4bool G4MCTSimVertex::GetStoreFlag() const
{
  return storeFlag;
}

inline void G4MCTSimVertex::SetInParticle(const G4MCTSimParticle* in)
{
  inParticleTrackID = in->GetTrackID();
}

inline void G4MCTSimVertex::SetInParticle(G4int in)
{
  inParticleTrackID = in;
}

inline G4int G4MCTSimVertex::GetInParticleTrackID() const
{
  return inParticleTrackID;
}

inline G4int G4MCTSimVertex::GetNofOutParticles() const
{
  return (G4int)outParticleTrackIDList.size();
}

inline G4int G4MCTSimVertex::AddOutParticle(const G4MCTSimParticle* out)
{
  outParticleTrackIDList.push_back(out->GetTrackID());
  return (G4int)outParticleTrackIDList.size();
}

inline G4int G4MCTSimVertex::AddOutParticle(G4int out)
{
  outParticleTrackIDList.push_back(out);
  return (G4int)outParticleTrackIDList.size();
}

inline G4int G4MCTSimVertex::GetOutParticleTrackID(G4int i) const
{
  G4int size = (G4int)outParticleTrackIDList.size();
  if(i >= 0 && i < size)
    return outParticleTrackIDList[i];
  else
    return 0;
}

#endif
