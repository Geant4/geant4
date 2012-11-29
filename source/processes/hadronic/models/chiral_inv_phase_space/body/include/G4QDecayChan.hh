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
//      ---------------- G4QDecayChan ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Decay Channel definition in CHIPS Model
// ------------------------------------------------------------
// Short description: In the CHIPS World the particles (G4QParticle)
// are defined. For unstable particles there is a G4QDecayChannelVector
// which describes different channels of decay for the particular
// particle. So the G4QDecayChannel class is the class for the description
// of such a decay channel in two or three particles (the secondaries can
// be unstable too and have firther decay).
// -------------------------------------------------------------------

#ifndef G4QDecayChan_h
#define G4QDecayChan_h 1

#include <iostream>
#include "globals.hh"
#include "G4QPDGCodeVector.hh"

class G4QDecayChan
{
public:
  // Constructors
  G4QDecayChan();                                                    // Default Constructor
  //G4QDecayChan(G4QPDGCodeVector secHadr, G4double probLimit = 1.); // General Constructor
  G4QDecayChan(G4double pLev, G4int PDG1, G4int PDG2, G4int PDG3=0); // DetailedConstructor
  G4QDecayChan(const G4QDecayChan& right);                       // Copy ByValueConstructor
  G4QDecayChan(G4QDecayChan* right);                           // Copy ByPointerConstructor

  ~G4QDecayChan();                                                   // Public Destructor

  // Operators
  const G4QDecayChan& operator=(const G4QDecayChan& right);
  G4bool               operator==(const G4QDecayChan& rhs) const;
  G4bool               operator!=(const G4QDecayChan& rhs) const;

  // Selectors
  G4double         GetDecayChanLimit() const;   // Get a Decay Channel Probability Limit
  G4double         GetMinMass() const;          // Get a Minimum Mass for the Decay Channel
  G4QPDGCodeVector GetVecOfSecHadrons();        // Get a Vector of secondary PDG-particles

  // Modifiers
  void SetDecayChanLimit(G4double newDecChanLim);// Set a Decay Channel Probability Limit
  void SetMinMass(G4double newMinMass);          // Set a Minimum Mass for the DecayChannel
  void SetVecOfSecHadrons(G4QPDGCodeVector hadV);// Set a Vector of secondary PDG-particles

  //private:
  // Encapsulated functions

private:
  // the Body
  G4double aDecayChanLimit;
  G4double theMinMass;
  G4QPDGCodeVector aVecOfSecHadrons;
};

// Not member operators
std::ostream&   operator<<(std::ostream& lhs, G4QDecayChan& rhs);
//----------------------------------------------------------------------------------------

inline G4bool G4QDecayChan::operator==(const G4QDecayChan& rhs) const {return this==&rhs;}
inline G4bool G4QDecayChan::operator!=(const G4QDecayChan& rhs) const {return this!=&rhs;}
 
inline G4double G4QDecayChan::GetDecayChanLimit() const    {return aDecayChanLimit;}
inline G4double G4QDecayChan::GetMinMass() const           {return theMinMass;}
inline G4QPDGCodeVector G4QDecayChan::GetVecOfSecHadrons() {return aVecOfSecHadrons;}

inline void G4QDecayChan::SetDecayChanLimit(G4double DCL)          {aDecayChanLimit=DCL;}
inline void G4QDecayChan::SetMinMass(G4double minMass)             {theMinMass=minMass;}
inline void G4QDecayChan::SetVecOfSecHadrons(G4QPDGCodeVector VSH) {aVecOfSecHadrons=VSH;}

#endif



