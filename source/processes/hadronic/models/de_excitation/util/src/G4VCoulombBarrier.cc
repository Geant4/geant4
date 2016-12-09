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
// $Id: G4VCoulombBarrier.cc 100690 2016-10-31 11:25:43Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)

#include "G4VCoulombBarrier.hh"
#include "G4PhysicalConstants.hh"

G4VCoulombBarrier::G4VCoulombBarrier(G4int anA, G4int aZ)
{
  theA = anA;
  theZ = aZ;
  theRho = 0.0; 
  theR0 = 1.5*CLHEP::fermi;
}

G4VCoulombBarrier::~G4VCoulombBarrier()
{}

void G4VCoulombBarrier::SetParameters(G4double rho, G4double r0)
{
  theRho = rho;
  theR0 = r0;
}
