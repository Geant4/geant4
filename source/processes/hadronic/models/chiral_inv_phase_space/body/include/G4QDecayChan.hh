// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QDecayChan.hh,v 1.2 2000-09-10 13:58:55 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QDecayChan_h
#define G4QDecayChan_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QDecayChan ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Decay Channel definition in CHIPS Model
// ------------------------------------------------------------

#include <iostream.h>
#include "globals.hh"
#include "G4QPDGCodeVector.hh"

class G4QDecayChan
{
public:
  // Constructors
  G4QDecayChan();                                                    // Default Constructor
  G4QDecayChan(G4QPDGCodeVector secHadr, G4double probLimit = 1.);   // General Constructor
  G4QDecayChan(G4double pLev, G4int PDG1, G4int PDG2, G4int PDG3=0); // Detailed Constructor
  G4QDecayChan(const G4QDecayChan& right);                           // Copy Constructor by value
  G4QDecayChan(G4QDecayChan* right);                                 // Copy Constructor by pointer

  ~G4QDecayChan();                                                   // Destructor

  // Operators
  const G4QDecayChan& operator=(const G4QDecayChan& right);
  G4int               operator==(const G4QDecayChan& rhs) const;
  G4int               operator!=(const G4QDecayChan& rhs) const;

  // Selectors
  G4double         GetDecayChanLimit() const;    // Get a Decay Channel Probability Limit
  G4double         GetMinMass() const;           // Get a Minimum Mass for the Decay Channel
  G4QPDGCodeVector GetVecOfSecHadrons();         // Get a Vector of secondary PDG-particles

  // Modifiers
  void SetDecayChanLimit(G4double newDecChanLim);// Set a Decay Channel Probability Limit
  void SetMinMass(G4double newMinMass);          // Set a Minimum Mass for the Decay Channel
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
ostream&   operator<<(ostream& lhs, G4QDecayChan& rhs);
//----------------------------------------------------------------------------------------

inline G4int G4QDecayChan::operator==(const G4QDecayChan& rhs) const {return this==&rhs;}
inline G4int G4QDecayChan::operator!=(const G4QDecayChan& rhs) const {return this!=&rhs;}
 
inline G4double G4QDecayChan::GetDecayChanLimit() const    {return aDecayChanLimit;}
inline G4double G4QDecayChan::GetMinMass() const           {return theMinMass;}
inline G4QPDGCodeVector G4QDecayChan::GetVecOfSecHadrons() {return aVecOfSecHadrons;}

inline void G4QDecayChan::SetDecayChanLimit(G4double DCL)          {aDecayChanLimit=DCL;}
inline void G4QDecayChan::SetMinMass(G4double minMass)             {theMinMass=minMass;}
inline void G4QDecayChan::SetVecOfSecHadrons(G4QPDGCodeVector VSH) {aVecOfSecHadrons=VSH;}

#endif



