// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmixt.cc,v 1.4 1999-05-18 02:40:46 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <iomanip.h>
#include <math.h>
#include "globals.hh"
#include "G3toG4.hh"
#include "G3EleTable.hh"
#include "G3MatTable.hh"
#include "G4Material.hh"

void PG4gsmixt(RWCString tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsmixt);

    // interpret the parameters
    G4String name = Spar[0].data();
    G4int imate = Ipar[0];
    G4int nlmat = Ipar[1];
    G4double dens = Rpar[0]*g/cm3;
    G4double *a = Rpar + 1;
    G4double *z = Rpar + 1+abs(nlmat);
    G4double *wmat = Rpar + 1 + 2*abs(nlmat);

    for (int i=0; i<abs(nlmat); i++){
      Rpar[i]=Rpar[i]*g/mole;
    };
    G4gsmixt(imate,name,a,z,dens,nlmat,wmat);
}

void G4gsmixt(G4int imate, G4String name, G4double a[], G4double z[],
              G4double dens, G4int nlmat, G4double wmat[]){
  G4int nmate = abs(nlmat);
  G4String sname = name.strip(RWCString::both);
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




