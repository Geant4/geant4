// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCHIPSWorld.hh,v 1.1 2000-09-04 07:46:20 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QCHIPSWorld_h
#define G4QCHIPSWorld_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCHIPSWorld ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for CHIPS World of particles in CHIPS Model
// ------------------------------------------------------------

#include <iostream.h>
#include "globals.hh"
#include "G4QParticleVector.hh"

class G4QCHIPSWorld
{
public:
  // Constructors
  G4QCHIPSWorld(G4int nOfParts = 0);                // Construction with N particles

  ~G4QCHIPSWorld();                                 // Destructor

  // Selectors
  G4QParticle* GetQParticle(G4int PDG)       const; // Get pointer to particle in CHIPS World
  G4QParticle* GetQParticle(G4QPDGCode QPDG) const; // Get pointer to particle in CHIPS World
  G4QParticle* GetQParticle(G4QPDGCode* pQP) const; // Get pointer to particle in CHIPS World
  G4int        GetQPEntries()                const; // Get a#of particles in CHIPS World

private:
  G4QParticleVector* InitCHIPSWorld(G4int nOfParts);// nOfParts<0 kills the CHIPS World @@??

private:
  G4QParticleVector* qWorld;                        // Pointer to the CHIPS World (C++ Hip)
};

 
inline G4QParticle* G4QCHIPSWorld::GetQParticle(G4int       PDG) const
{
  return (*qWorld)[G4QPDGCode(PDG).GetQCode()];
}
inline G4QParticle* G4QCHIPSWorld::GetQParticle(G4QPDGCode QPDG) const
{
  return (*qWorld)[QPDG.GetQCode()];
}
inline G4QParticle* G4QCHIPSWorld::GetQParticle(G4QPDGCode* pQP) const
{
  return (*qWorld)[pQP->GetQCode()];
}
inline G4int        G4QCHIPSWorld::GetQPEntries() const {return (*qWorld).entries();}

#endif



