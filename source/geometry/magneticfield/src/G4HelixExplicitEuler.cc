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
// $Id: G4HelixExplicitEuler.cc,v 1.6 2006/06/29 18:24:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//
//  Helix Explicit Euler: x_1 = x_0 + helix(h)
//  with helix(h) being a helix piece of length h
//  most simple approach for solving linear differential equations.
//  Take the current derivative and add it to the current position.
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// -------------------------------------------------------------------

#include "G4HelixExplicitEuler.hh"
#include "G4ThreeVector.hh"

void
G4HelixExplicitEuler::DumbStepper( const G4double  yIn[],
				   G4ThreeVector   Bfld,
				   G4double  h,
				   G4double  yOut[])
{
  AdvanceHelix(yIn, Bfld, h, yOut);

  // NormaliseTangentVector( yOut );  // this could harm more than it helps 
}  
