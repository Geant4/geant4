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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4QCHIPSWorld.hh,v 1.10 2001-11-26 14:11:45 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCHIPSWorld ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for CHIPS World of particles in CHIPS Model
// ------------------------------------------------------------

#ifndef G4QCHIPSWorld_h
#define G4QCHIPSWorld_h 1

#include "g4std/iostream"
#include "globals.hh"
#include "G4QParticleVector.hh"

class G4QCHIPSWorld
{
public:
  // Constructors
  G4QCHIPSWorld(G4int nOfParts = 0);                // Construction with N particles
  G4QCHIPSWorld(const G4QCHIPSWorld& right);        // Copy Constructor by value
  G4QCHIPSWorld(G4QCHIPSWorld* right);              // Copy Constructor by pointer

  ~G4QCHIPSWorld();                                 // Destructor

  // Overloaded Operators
  const G4QCHIPSWorld& operator=(const G4QCHIPSWorld& right);

  // Selectors
  G4QParticle* GetQParticle(G4int PDG)       const; // Get pointer to particle in CHIPS World
  G4QParticle* GetQParticle(G4QPDGCode QPDG) const; // Get pointer to particle in CHIPS World
  G4QParticle* GetQParticle(G4QPDGCode* pQP) const; // Get pointer to particle in CHIPS World
  G4int        GetQPEntries()                const; // Get a#of particles in CHIPS World

private:
  G4QParticleVector* InitCHIPSWorld(G4int nOfParts);// nOfParts<0 kills the CHIPS World @@??

// Body
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
inline G4int        G4QCHIPSWorld::GetQPEntries() const {return (*qWorld).size();}

#endif



