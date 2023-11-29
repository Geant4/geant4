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
//
// 
// John Allison  30th October 1996

// Class Description:
// Base class for all things visible, i.e., which have Vis Attributes.
//
// Note: a null pointer implies no attributes.  Under those circumstances
// the visualization system is free to choose some.
// Class Description - End:


#ifndef G4VISIBLE_HH
#define G4VISIBLE_HH

#include "globals.hh"
#include <iostream>

class G4VisAttributes;

class G4Visible {

  friend std::ostream& operator << (std::ostream& os, const G4Visible& v);

public: // With description

  G4Visible ();
  G4Visible (const G4Visible&);
  G4Visible (G4Visible&&);
  G4Visible (const G4VisAttributes*);

  virtual ~G4Visible ();

  G4Visible& operator= (const G4Visible&);
  G4Visible& operator= (G4Visible&&);

  G4bool operator != (const G4Visible& right) const;

  const G4VisAttributes* GetVisAttributes () const;

  void SetVisAttributes (const G4VisAttributes*);
  // The G4VisAttributes object is not stored in a G4Visible; only a
  // reference, a const pointer, is kept.  Therefore the
  // G4VisAttributes object to which it refers must have a life long
  // enough to satisfy all uses of the G4Visible object.  E.g., if the
  // G4Visible object is created on the heap (using `new') then the
  // associated G4VisAttributes object would normally also be created
  // on the heap and managed in the same way.

  void SetVisAttributes (const G4VisAttributes&);
  // A copy of the G4VisAttributes object is created on the heap to
  // ensure a long life.

  // Access functions to the string for user customizable information
  virtual   const G4String&  GetInfo() const;
  virtual   void             SetInfo(const G4String& info);

protected:

  G4String  fInfo;   // String for user customizable information
  const G4VisAttributes* fpVisAttributes;
  G4bool fAllocatedVisAttributes;
};

#include "G4Visible.icc"

#endif
