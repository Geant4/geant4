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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)
//
// 14-11-2007 modified barrier by JMQ (test30) 
// 15-11-2010 V.Ivanchenko use G4Pow and cleanup 

#include "G4FermiCoulombBarrier.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4NuclearRadii.hh"

G4FermiCoulombBarrier::G4FermiCoulombBarrier(G4int A, G4int Z)
  : G4VCoulombBarrier(A, Z)
{
  SetParameters(G4NuclearRadii::RadiusCB(Z, A), 1.3*CLHEP::fermi);
  factor = CLHEP::elm_coupling*0.6*g4calc->Z13(7)/theRho;
}

G4double G4FermiCoulombBarrier::GetCoulombBarrier(
         G4int ARes, G4int ZRes, G4double) const 
{
  if(0 == theZ) { return 0.0; }
  G4int A = theA + ARes;
  G4int Z = theZ + ZRes;
  G4double cb = factor*((Z*Z)/g4calc->Z13(A) - (theZ*theZ)/g4calc->Z13(theA)
			- (ZRes*ZRes)/g4calc->Z13(ARes));
  return cb;
}

