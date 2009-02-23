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
// $Id: G4QParticle.hh,v 1.25 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QParticle ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Particles in the CHIPS Model
// ---------------------------------------------------
// Short description: The G4QParticle is a part of the CHIPS World. It is
// characterized by the quark content, spin, mass, width and a vector of
// the decay channels (G4QDecayCannelVector).
// -----------------------------------------------------------------------

#ifndef G4QParticle_h
#define G4QParticle_h 1

#include <iostream>
#include "globals.hh"
#include "G4QDecayChanVector.hh"

class G4QParticle
{
public:
  // Constructors
  G4QParticle();                             // Default Constructor
  G4QParticle(G4bool f, G4int theQCode);     // QCode Constructor, f-verbose
  G4QParticle(G4int thePDG);                 // PDGCode Constructor
  G4QParticle(const G4QParticle& right);     // Copy Constructor by value
  G4QParticle(G4QParticle* right);           // Copy Constructor by pointer

  ~G4QParticle();                            // Public Destructor

  // Operators
  const G4QParticle& operator=(const G4QParticle& right);
  G4bool             operator==(const G4QParticle& rhs) const;
  G4bool             operator!=(const G4QParticle& rhs) const;

  // Selectors
  G4QPDGCode          GetQPDG() const;          // Get a PDG-Particle of the Particle
  G4int               GetPDGCode() const;       // Get a PDG Code of the Particle
  G4int               GetQCode() const;         // Get a Q Code of the Particle
  G4int               GetSpin() const;          // Get 2s+1 of the Particle
  G4int               GetCharge() const;        // Get a Charge of the Particle
  G4int               GetStrange() const;       // Get a Strangeness of the Particle
  G4int               GetBaryNum() const;       // Get a Baryon Number of the Particle
  G4QContent          GetQContent();            // Get Quark Content of the Particle
  G4QDecayChanVector  GetDecayVector();         // Get a Decay Vector for the Particle
  G4double            GetMass();                // Get a mass value for the Particle
  G4double            GetWidth();               // Get a width value for the Particle

  // Modifiers
  G4QDecayChanVector InitDecayVector(G4int Q);// Init DecayVector in theCHIPSWorld by QCode
  void InitPDGParticle(G4int thePDGCode);
  void InitQParticle(G4int theQCode);

  // General
  G4double MinMassOfFragm();                    // Minimal mass of decaing fragments

private:
  // Encapsulated functions

private:
  // the Body
  G4QPDGCode          aQPDG;
  G4QDecayChanVector  aDecay;
  G4QContent          aQuarkCont;       // @@ Secondary (added for acceleration - check)
};

// Not member operators
std::ostream&   operator<<(std::ostream& lhs, G4QParticle& rhs);
// Not member functions
//----------------------------------------------------------------------------------------

inline G4bool G4QParticle::operator==(const G4QParticle& rhs) const {return this==&rhs;}
inline G4bool G4QParticle::operator!=(const G4QParticle& rhs) const {return this!=&rhs;}
 
inline G4QPDGCode    G4QParticle::GetQPDG()    const {return aQPDG;}
inline G4int         G4QParticle::GetQCode()   const {return aQPDG.GetQCode();}
inline G4int         G4QParticle::GetPDGCode() const {return aQPDG.GetPDGCode();}
inline G4int         G4QParticle::GetSpin()    const {return aQPDG.GetSpin();}
inline G4int         G4QParticle::GetCharge()  const {return aQuarkCont.GetCharge();}
inline G4int         G4QParticle::GetStrange() const {return aQuarkCont.GetStrangeness();}
inline G4int         G4QParticle::GetBaryNum() const {return aQuarkCont.GetBaryonNumber();}
inline G4QContent    G4QParticle::GetQContent()      {return aQuarkCont;}
inline G4QDecayChanVector G4QParticle::GetDecayVector() {return aDecay;}
inline G4double      G4QParticle::GetMass()          {return aQPDG.GetMass();}
inline G4double      G4QParticle::GetWidth()         {return aQPDG.GetWidth();}

inline G4double G4QParticle::MinMassOfFragm()
{
  G4int nCh=aDecay.size();
  G4double m=GetMass();
  G4double min=m;
  if(nCh)
  {
    min=aDecay[0]->GetMinMass();
    if(nCh>1) for(G4int j=1; j<nCh; j++)
    {
      G4double next=aDecay[j]->GetMinMass();
      if(next<min) min=next;
    }
  }
  G4double w=GetWidth();
  G4double lim=m+.001;
  if(w)   lim-=1.5*w;
  if(min<lim) min=lim;
  return min;
}

#endif



