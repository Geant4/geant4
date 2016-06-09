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
// $Id: G4QBesIKJY.hh,v 1.3 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QBesIKJY ----------------
//             by Mikhail Kossov, Sept 2000.
//  class header for Bessel I0/I1 and K0/K1 functions in CHIPS Model
// --------------------------------------------------------------------
// Short description: CHIPS bessel functions for mass and scattering
// integrals.
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
