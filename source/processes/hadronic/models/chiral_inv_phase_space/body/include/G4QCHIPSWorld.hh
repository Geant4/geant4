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
//
// $Id$
//
//      ---------------- G4QCHIPSWorld ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for CHIPS World of particles in CHIPS Model
// ------------------------------------------------------------
// Short description: The CHIPS World is a world of elementary particles
// and nuclear fragments. This class is a singletone, but without fixed
// limits. E.g. the nuclear fragments as possible G4Candidates can be
// initialized in the CHIPS World only up to Be8 od C12 or other bigger
// nuclear fragment. If one need the heavy fragment production then the
// the CHIPS World must be initialized up to these fragments (see the
// CHIPS Manual), but the price in performans will be big, because in
// each act of the fragmentation competition these numerous candidates
// take place in the competition and the hadronization probability is
// calculated each time for each of them, so the Be8 limit (Be8->alpha+
// alpha decays very fast and contribute to the alpha-spectrum) is the
// most optimal.
// -------------------------------------------------------------------

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



