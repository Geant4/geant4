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

#include "G4NuclearFermiDensity.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"

G4NuclearFermiDensity::G4NuclearFermiDensity(G4int anA, G4int /*aZ*/) 
 : theA(anA), a(0.545 * fermi) 
{
  G4double a13 = G4Pow::GetInstance()->Z13(anA);
  const G4double r0 = 1.16 * (1. - 1.16/(a13*a13)) * fermi;
  theR = r0 * a13;
  Setrho0(3./ (4.*pi *r0*r0*r0 * theA * (1. + sqr(a/theR)*pi2 )));
}

G4NuclearFermiDensity::~G4NuclearFermiDensity() {}

