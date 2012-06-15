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
// $Id: G4RandomDirection.hh,v 1.5 2008-03-19 17:00:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
// Class description:
//
// Function returning a unit 3-vector homogeneously randomised over 4pi
// solid angle. It can be used in any particle scattering methods
// instead of:
//   z=R1, x=SQRT(1-R1*R1)*SIN(2*pi*R2), y=SQRT(1-R1*R1)*COS(2*pi*R2)
// providing more performant results.

// History:
//    18.03.08 V. Grichine, unit radius sphere surface based algorithm
//      ~ 2007 M. Kossov, algorithm based on 8 Quadrants technique
//
// ------------------------------------------------------------
#ifndef G4RANDOMDIR_HH
#define G4RANDOMDIR_HH

#include <CLHEP/Units/PhysicalConstants.h>
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

inline G4ThreeVector G4RandomDirection()
{
  G4double cosTheta  = 2.*G4UniformRand()-1.;
  G4double sinTheta2 = 1. - cosTheta*cosTheta;
  if( sinTheta2 < 0.)  sinTheta2 = 0.;
  G4double sinTheta  = std::sqrt(sinTheta2); 
  G4double phi       = CLHEP::twopi*G4UniformRand();
  return G4ThreeVector(sinTheta*std::cos(phi),
                       sinTheta*std::sin(phi), cosTheta).unit(); 
}

#endif  /* G4RANDOMDIR_HH */
