// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmate.cc,v 1.1 1999-01-07 16:06:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------------
// History:
// 25-Feb-1997 Lockman - removed isMixture when calling G4Material constr.
//                       to conform with P. Maire's new interface 29-Jan-1997
//                     - added units to density, atomic weight
//                     - added check for Vacuum
// 14-Jun-1998 Lockman - all states initialized to kStateUndefined;
//                       logic to handle G4 vacuum state modified.
// ---------------------------------------------------------------------------

#include <math.h>
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"

//extern int debugOn;

void PG4gsmate(RWCString tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsmate);
    G4String name = Spar[0];
    G4int imate = Ipar[0];
    G4int nwbf = Ipar[1];
    G4double a = Rpar[0];
    G4double z = Rpar[1];
    G4double dens = Rpar[2];
    G4double radl = Rpar[3];
    G4double absl = Rpar[4];
    G4double *ubuf = &Rpar[5];

    G4gsmate(imate,name,a,z,dens,radl,nwbf,ubuf);
}

void G4gsmate(G4int imate, G4String name, G4double ain, G4double zin,
              G4double densin, G4double radl, G4int nwbf, G4double* ubuf)
{

        // set default arguments
    
    G4double zdef = 1.;
    G4double adefval = 2.;
    G4double adef = adefval*g/mole;
    G4double G3_minimum_density = 1.e-10*g/cm3;
    G4double a = ain*g/mole;
    G4double z = zin;
    G4double dens = densin*g/cm3;
    
    if (dens < G3_minimum_density) {
      dens=G3_minimum_density;
      G4cerr << "G3 vacuum state encountered for material " << name 
	   << ". Set dens = " << dens*cm3/g << " g/cm3" << endl;
    }
    if (z<1) {
      z=zdef;
      a=adef;
      G4cerr << "Z<1 (" << z << ") for material " << name 
	   << ". Set Z = " << z << " A = " << adefval << " g/mole" << endl;
    }
    G4Material* mte = new G4Material(name, z, a, dens);
    
        // add the material to the List
    G3Mat.put(&imate,mte);
}







