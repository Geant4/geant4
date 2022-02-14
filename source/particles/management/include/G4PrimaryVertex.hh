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
// G4PrimaryVertex
//
// class description:
//
// This is the class which represents a primary vertex. The object of this
// class is set to G4Event by a G4VPrimaryGenerator concrete class.
// This class object has one or more G4PrimaryParticle objects as primary
// particles.

// Authors: G.Cosmo, 2 December 1995 - Design, based on object model
//          M.Asai, 29 January 1996 - First implementation
// --------------------------------------------------------------------
#ifndef G4PrimaryVertex_hh
#define G4PrimaryVertex_hh 1

#include "globals.hh"
#include "pwdefs.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4PrimaryParticle.hh"

class G4VUserPrimaryVertexInformation;

class G4PrimaryVertex 
{
  public:

    G4PrimaryVertex();
    G4PrimaryVertex(G4double x0,G4double y0,G4double z0,G4double t0);
    G4PrimaryVertex(G4ThreeVector xyz0, G4double t0);
      // Constructors

    virtual ~G4PrimaryVertex();
      // Destructor

    G4PrimaryVertex(const G4PrimaryVertex& right);
    G4PrimaryVertex& operator=(const G4PrimaryVertex& right);
      // Copy constructor & assignment operator

    G4bool operator==(const G4PrimaryVertex& right) const;
    G4bool operator!=(const G4PrimaryVertex& right) const;
      // Equality operators

    inline void* operator new(size_t);
    inline void operator delete(void* aPrimaryVertex);
      // Overloaded new/delete operators

    // Accessors & modifiers

    inline G4ThreeVector GetPosition() const;
    inline void SetPosition(G4double x0, G4double y0, G4double z0);
    inline G4double GetX0() const;
    inline G4double GetY0() const;
    inline G4double GetZ0() const;
    inline G4double GetT0() const;
    inline void SetT0(G4double t0);
    inline G4int GetNumberOfParticle() const;
    inline void SetPrimary(G4PrimaryParticle* pp);
    G4PrimaryParticle* GetPrimary(G4int i=0) const;
    inline void SetNext(G4PrimaryVertex* nv);
    inline void ClearNext();
    inline G4PrimaryVertex* GetNext() const;
    inline G4double GetWeight() const;
    inline void SetWeight(G4double w);
    inline void SetUserInformation(G4VUserPrimaryVertexInformation* info);
    inline G4VUserPrimaryVertexInformation* GetUserInformation() const;

    void Print() const;
  
 private:

    G4double X0 = 0.0;
    G4double Y0 = 0.0;
    G4double Z0 = 0.0;
    G4double T0 = 0.0;
    G4PrimaryParticle* theParticle = nullptr;
    G4PrimaryParticle* theTail = nullptr;
    G4PrimaryVertex* nextVertex = nullptr;
    G4PrimaryVertex* tailVertex = nullptr;
    G4double Weight0 = 1.0;
    G4VUserPrimaryVertexInformation* userInfo = nullptr;
    G4int numberOfParticle = 0;
};

extern G4PART_DLL G4Allocator<G4PrimaryVertex>*& aPrimaryVertexAllocator();

// ------------------------
// Inline methods
// ------------------------

inline
void* G4PrimaryVertex::operator new(std::size_t)
{
  if (aPrimaryVertexAllocator() == nullptr)
  {
    aPrimaryVertexAllocator() = new G4Allocator<G4PrimaryVertex>;
  }
  return (void*) aPrimaryVertexAllocator()->MallocSingle();
}

inline
void G4PrimaryVertex::operator delete(void* aPrimaryVertex)
{
  aPrimaryVertexAllocator()->FreeSingle((G4PrimaryVertex*) aPrimaryVertex);
}

inline
G4ThreeVector  G4PrimaryVertex::GetPosition() const
{
  return G4ThreeVector(X0,Y0,Z0);
}

inline
void G4PrimaryVertex::SetPosition(G4double x0,G4double y0,G4double z0)
{
  X0 = x0; Y0 = y0; Z0 = z0;
}

inline
G4double G4PrimaryVertex::GetX0() const
{
  return X0;
}

inline
G4double G4PrimaryVertex::GetY0() const
{
  return Y0;
}

inline
G4double G4PrimaryVertex::GetZ0() const
{
  return Z0;
}

inline
G4double G4PrimaryVertex::GetT0() const
{
  return T0;
}

inline
void G4PrimaryVertex::SetT0(G4double t0)
{
  T0 = t0;
}

inline
G4int G4PrimaryVertex::GetNumberOfParticle() const
{
  return numberOfParticle;
}

inline
void G4PrimaryVertex::SetPrimary(G4PrimaryParticle* pp)
{ 
  if(theParticle == nullptr) { theParticle = pp;     }
  else                       { theTail->SetNext(pp); }
  theTail = pp;
  ++numberOfParticle;
}

inline
void G4PrimaryVertex::SetNext(G4PrimaryVertex* nv)
{ 
  if(nextVertex == nullptr) { nextVertex = nv; }
  else                      { tailVertex->SetNext(nv); }
  tailVertex = nv;
}

inline
void G4PrimaryVertex::ClearNext()
{
  nextVertex = nullptr;
  tailVertex = nullptr;
}

inline
G4PrimaryVertex* G4PrimaryVertex::GetNext() const
{
  return nextVertex;
}

inline
G4double G4PrimaryVertex::GetWeight() const
{
  return Weight0;
}

inline
void G4PrimaryVertex::SetWeight(G4double w)
{
  Weight0 = w;
}

inline
void G4PrimaryVertex::SetUserInformation(G4VUserPrimaryVertexInformation* info)
{
  userInfo = info;
}

inline
G4VUserPrimaryVertexInformation* G4PrimaryVertex::GetUserInformation() const
{
  return userInfo;
}

#endif

