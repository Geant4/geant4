// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmixt.cc,v 1.3 1999-05-06 17:46:59 lockman Exp $
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
  for (int i=0; i< nmate; i++){
    G4Element* theElement = G3Ele.GetEle(z[i]);
    if (nlmat>0) {
      G4double fractionmass = wmat[i];
      theMixture->AddElement(theElement, fractionmass);
    } else if (nlmat<0) {
      G4int natoms = wmat[i];
      theMixture->AddElement(theElement, natoms);
    }
  }
  G3Mat.put(imate, theMixture);
}




