// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmate.cc,v 1.3 1999-11-15 10:39:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G3toG4.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G3MatTable.hh"

//extern int debugOn;

void PG4gsmate(G4String tokens[])
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

  G4gsmate(imate, name, a, z, dens, radl, nwbf, ubuf);
}

void G4gsmate(G4int imate, G4String name, G4double ain, G4double zin,
              G4double densin, G4double radl, G4int nwbf, G4double* ubuf)
{

  // set default arguments
    
  G4double zdef = 1.;
  G4double adef = 1.01*g/mole;
  G4double G3_minimum_density = 1.e-10*g/cm3;
  G4double theA = ain*g/mole;
  G4double theZ = zin;
  G4double theDensity = densin*g/cm3;

  G4Material* mte = 0;
  G4State theState = kStateUndefined;
  G4double thePressure = STP_Pressure;
  G4double theTemperature = STP_Temperature;
  G4String sname = name.strip(G4String::both);
  G4String symbol;

  if (sname == "AIR") {

    // handle the built in AIR mixture
    G4double a[2], z[2], wmat[2];
    a[0] = 14.01*g/mole;
    a[1] = 16.00*g/mole;
    z[0] = 7;
    z[1] = 8;
    wmat[0] = 0.7;
    wmat[1] = 0.3;
    theDensity = 1.2931*mg/cm3;
    int n=2;
    G4gsmixt(imate, sname, a, z, theDensity, n, wmat);
  } else {
    if (theDensity < G3_minimum_density) {

    // handle the built in VACUUM
    theDensity = universe_mean_density;
    theState = kStateGas;
    theZ = zdef;
    theA = adef;
    thePressure = 3.e-18*pascal;
    theTemperature = 2.73*kelvin;

    G4cout << "G3 vacuum " 
	   << sname << " (Z " << zin << ", A " << ain << " g/mole, density " 
	   << densin << " g/cm3)" << endl
	   << "Creating G4 vacuum "
	   << "(Z " << theZ << ", A " << theA*mole/g << ", g/mole, "
	   << " density " << theDensity*cm3/g << " g/cm3, "
	   << " pressure " << thePressure/pascal
	   << " pascal, temp " << theTemperature/kelvin 
	   << " kelvins)" << endl;
    }
    G4Material* mte = new 
      G4Material(sname, theZ, theA, theDensity, theState, theTemperature, 
		 thePressure);
    G3Mat.put(imate, mte);
  }
}







