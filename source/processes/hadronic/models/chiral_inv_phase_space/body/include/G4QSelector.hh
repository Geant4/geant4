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
// $Id: G4QSelector.hh,v 1.7 2001-11-26 14:11:45 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QSelector ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for decay Selection for Hadrons of CHIPS Model
// ------------------------------------------------------------

#ifndef G4QSelector_h
#define G4QSelector_h 1

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
