// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QSelector.hh,v 1.1 2000-08-17 13:55:08 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QSelector_h
#define G4QSelector_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QSelector ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for decay Selection for Hadrons of CHIPS Model
// ------------------------------------------------------------

//#include "globals.hh"
#include "G4QHadronVector.hh"
#include "G4QPDGCodeVector.hh"

class G4QSelector
{
public:
  // Constructors
  G4QSelector();                                      // Default Constructor
  G4QSelector(G4QHadronVector hadronVec);             // Constructor for the Hadrons Set
  ~G4QSelector();                                     // Destructor
  // Operators
  const G4QSelector& operator=(const G4QSelector& right);
  G4int operator==(const G4QSelector& right) const;
  G4int operator!=(const G4QSelector& right) const;
  // Selectors
  G4QHadronVector GetHadrons()    const;              // Get PDG code of the Hadron
  // Modifiers
  void   SetHadrons(G4QHadronVector hadronVec);       // Get PDG code of the Hadron
  G4bool SelectPDGSet(G4QPDGCodeVector thePDGCodes);  // Select event with the PDGCode Set
private:  
  G4QHadronVector theHadrons;                         // Vector of Hadron
};

inline G4int G4QSelector::operator==(const G4QSelector& right) const {return this==&right;}
inline G4int G4QSelector::operator!=(const G4QSelector& right) const {return this!=&right;}
 
inline G4QHadronVector G4QSelector::GetHadrons()             const   {return theHadrons;}

inline void G4QSelector::SetHadrons(G4QHadronVector hadronVec)       {theHadrons=hadronVec;}

#endif
