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
//

#include "G4NuclearFermiDensity.hh"

G4NuclearFermiDensity::G4NuclearFermiDensity(G4double anA, G4double aZ) 
  :  a(0.545 * fermi) 
{
//        const G4double r0=1.14*fermi;
	const G4double r0=1.16 * ( 1. - 1.16 * std::pow(anA, -2./3.)) * fermi;
	theA = G4int(anA);
	theZ = G4int(aZ);
	theR= r0 * std::pow(anA, 1./3. );
	Setrho0(3./ (4. * pi * std::pow(r0,3.) * theA * ( 1. + sqr(a/theR)*pi2 )));
}

G4NuclearFermiDensity::~G4NuclearFermiDensity() {}
        

