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
// $Id: G4QBesIKJY.hh,v 1.1 2006-06-02 13:38:27 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QBesIKJY ----------------
//             by Mikhail Kossov, Sept 2000.
//  class header for Bessel I0/I1 and K0/K1 functions in CHIPS Model
// --------------------------------------------------------------------

#ifndef G4QBesIKJY_h
#define G4QBesIKJY_h 1

#include <iostream>
#include "globals.hh"

// Ennumerated type
enum G4QBIType {BessI0, BessI1, EBessI0, EBessI1, BessJ0, BessJ1,
                BessK0, BessK1, EBessK0, EBessK1, BessY0, BessY1};
class G4QBesIKJY
{
public:
  // Constructor/Destructor
  G4QBesIKJY(G4QBIType type = BessI0); // the Default Construction is for BesselI0
  ~G4QBesIKJY();

public:
  // Member Functions
  G4double operator ()(G4double x) const;

private:
  // Body
  G4bool ex;
  G4bool ij;
  G4bool ik;
  G4int  nu;
};

#endif
