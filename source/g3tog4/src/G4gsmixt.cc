// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmixt.cc,v 1.6 1999-12-05 17:50:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "g4std/iomanip"
#include "g4std/strstream"
#include "g4std/iomanip"
#include <math.h>

#include "globals.hh"
#include "G3toG4.hh"
#include "G3EleTable.hh"
#include "G3MatTable.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"

void PG4gsmixt(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsmixt);

    // interpret the parameters
    G4String name = Spar[0].data();
    G4int imate = Ipar[0];
    G4int nlmat = Ipar[1];
    //G4double dens = Rpar[0]*g/cm3;
    G4double dens = Rpar[0];
    G4double *a = Rpar + 1;
    G4double *z = Rpar + 1+abs(nlmat);
    G4double *wmat = Rpar + 1 + 2*abs(nlmat);

    for (int i=0; i<abs(nlmat); i++){
      //Rpar[i]=Rpar[i]*g/mole;
      Rpar[i]=Rpar[i];
    };
    G4gsmixt(imate,name,a,z,dens,nlmat,wmat);
}

// replaced with G3EleTable 
// only used G4Elements are created;
// !! no checking of given A of the element;
//
// extern G4Element* CreateElement(G4double zeff, G4double aeff, G4String matName);


void G4gsmixt(G4int imate, G4String name, G4double* a, G4double* z,
              G4double dens, G4int nlmat, G4double* wmat)
{
  // in Geant3:
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  G4int i=0;
  if (nlmat<0) {
    // in case of proportions given in atom counts (nlmat<0),
    // the wmat[i] are converted to weight fractions
    G4double aMol = 0.;
    for (i=0; i<abs(nlmat); i++) { 
      // total molecular weight 
      aMol += wmat[i]*a[i]; 
    }  
    if (aMol == 0.)
      G4Exception("\nG4mixt: Total molecular weight in " +
        name + " = 0.");       
    for (i=0; i<abs(nlmat); i++) {
      // weight fractions
      wmat[i] = wmat[i]*a[i]/aMol;
    }
  }

  // create material with given number of components
  // (elements)

  G4Material* material 
    = new G4Material(name, dens*g/cm3, abs(nlmat));
  for (i=0; i<abs(nlmat); i++) {
    // add units
    // G4Element* element = G4Element(z[i], a[i]*g/mole, name);
    G4Element* element = G3Ele.GetEle(z[i]);
    material->AddElement(element, wmat[i]);    
  }

  // add the material to the List
  G3Mat.put(imate, material);
}

/*
void G4gsmixt(G4int imate, G4String name, G4double a[], G4double z[],
              G4double dens, G4int nlmat, G4double wmat[]){
  G4int nmate = abs(nlmat);
  G4String sname = name.strip(G4String::both);
  G4double theDensity = dens*g/cm3;

  G4Material* theMixture = new G4Material(name, dens, nmate); 
  G4bool ok=true;
  for (int i=0; i< nmate; i++){
    G4Element* theElement = G3Ele.GetEle(z[i]);
    if (nlmat>0) {
      G4double fractionmass = wmat[i];
      ok = ok && abs(fractionmass)<=1.;
      theMixture->AddElement(theElement, fractionmass);
    } else if (nlmat<0) {
      G4int natoms = wmat[i];
      ok = ok && wmat[i] == natoms;
      theMixture->AddElement(theElement, natoms);
    } else {
      ok=false;
    }
  }
  if (ok) {
    G3Mat.put(imate, theMixture);
  } else {
    if (nlmat>0) {
      G4cerr << "G4gsmixt: for mixture '" << name 
	     << "' some |weights|>1 : " << endl;
      for (G4int i=0;i<nlmat; i++) {
	G4cerr << "Component " << setw(3) << i+1 << " fraction: "
	       << setw(10) << wmat[i] << endl;
      }
    } else if (nlmat<0) {
      G4cerr << "G4gsmixt: for mixture '" << name 
	     << "' some #natoms are non-integer: " << endl;
      for (G4int i=0;i<nlmat; i++) {
	G4cerr << "Component " << setw(3) << i+1 << " #atoms "
	       << setw(10) << wmat[i] << endl;
      }
    } else {
      G4cerr << "G4gsmixt: Number of components for mixture '" 
	     << name << "' (" << nlmat << ") not allowed." << endl;
    }
  }
}
*/



