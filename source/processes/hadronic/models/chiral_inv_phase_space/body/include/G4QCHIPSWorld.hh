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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4QCHIPSWorld.hh,v 1.23 2004/12/01 15:16:28 mkossov Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//      ---------------- G4QCHIPSWorld ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for CHIPS World of particles in CHIPS Model
// ------------------------------------------------------------

#ifndef G4QCHIPSWorld_h
#define G4QCHIPSWorld_h 1

#include <iostream>
#include "globals.hh"
#include "G4QParticleVector.hh"

class G4QCHIPSWorld
{
  // Constructor/Destructor
protected:
  G4QCHIPSWorld();             // the Default Construction is protected - Singelton
public:
  ~G4QCHIPSWorld();            // Destructor is public because of Windows compilation error

  // Member Functions
private:
  G4QCHIPSWorld(const G4QCHIPSWorld& right);       // copy QCHIPSWorld by value
  G4QCHIPSWorld(G4QCHIPSWorld* right);             // copy QCHIPSWorld by pointer
  const G4QCHIPSWorld& operator=(const G4QCHIPSWorld& right);//copy QCHIPSWorld by Operator
  G4bool operator==(const G4QCHIPSWorld &right) const;
  G4bool operator!=(const G4QCHIPSWorld &right) const;

public:
  // Pointers to Particles of the Singeltone of the CHIPS World
  static G4QCHIPSWorld* Get();
  G4QParticleVector* GetParticles(G4int nOfParts=0);
  // Selectors
  G4QParticle* GetQParticle(G4int PDG)       const;// Get pointer to particle in CHIPSWorld
  G4QParticle* GetQParticle(G4QPDGCode QPDG) const;// Get pointer to particle in CHIPSWorld
  G4QParticle* GetQParticle(G4QPDGCode* pQP) const;// Get pointer to particle in CHIPSWorld
  G4int        GetQPEntries()                const;// Get a#of particles in CHIPS World

// Body
private:
  //static G4QCHIPSWorld* aWorld;             // Pointer to the CHIPS World
  static G4QParticleVector& GetQWorld();
};

 
inline G4QParticle* G4QCHIPSWorld::GetQParticle(G4int       PDG) const
{
  G4int qCode=G4QPDGCode(PDG).GetQCode();
  //G4cout<<"G4QCHIPSWorld::GetQPart:Q="<<qCode<<",Max="<<qWorld.size()<<G4endl;
  return GetQWorld()[qCode];
}

inline G4QParticle* G4QCHIPSWorld::GetQParticle(G4QPDGCode QPDG) const
{
  return GetQWorld()[QPDG.GetQCode()];
}

inline G4QParticle* G4QCHIPSWorld::GetQParticle(G4QPDGCode* pQP) const
{
  return GetQWorld()[pQP->GetQCode()];
}

inline G4int        G4QCHIPSWorld::GetQPEntries() const {return GetQWorld().size();}

#endif



