// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmixt.cc,v 1.1 1999-01-07 16:06:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------------
// History:
// 25-Feb-1997 Lockman - added units to density, atomic weight
// ---------------------------------------------------------------------------

#include "G4ios.hh"
#include <strstream.h>
#include <iomanip.h>
#include <math.h>
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
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
    for (int i=0; i<abs(nlmat); i++){
        Rpar[i]=Rpar[i]*g/mole;
    };
    G4double *a = &Rpar[1];
    G4double *z = &Rpar[1+abs(nlmat)];
    G4double *wmat = &Rpar[1+2*abs(nlmat)];

    G4gsmixt(imate,name,a,z,dens,nlmat,wmat);
}

void G4gsmixt(G4int imate, G4String name, G4double a[], G4double z[],
              G4double dens, G4int nlmat, G4double wmat[])
{
  // Build material components as 'elements'
  char indx[5], symbol[20];
  G4Element* el;
  G4String elName, elSymbol;
  G4bool isMixture = true;
  G4Material* gmate = new G4Material(name, dens, nlmat);

  G4double wt=0., zeff=0., aeff=0., frac;

      // for case of proportions given in atom counts (nlmat<0),
      // convert to weight fractions

  ostrstream ostr_indx(indx, sizeof indx);
  ostrstream ostr_symb(symbol, sizeof symbol);
  
  for (G4int i=0; i<abs(nlmat); i++) {
//    printf(indx,"%d\n",i);
//    printf(symbol,"Z%dA%d\n",int(z[i]),int(a[i]));
      ostr_indx << setw(3) << i;
      ostr_symb << setw(3) << int(z[i]) << " " << setw(3) << int(a[i]);
      elName = "Material "+name+" component "+indx;
      elSymbol = symbol;
      G4cout << "elName: " << elName << endl;
      G4cout << "elSymbol: " << elSymbol << endl;

          //
          // mod 25-Feb-1997 Lockman, removed isMixture to match
          // G4Element interface
          //
      
//    el = new G4Element(elName, elSymbol, z[i], a[i], isMixture);
    el = new G4Element(elName, elSymbol, z[i], a[i]);
    if ( nlmat < 0 ) {
      gmate->AddElement(el, int(wmat[i]));
    } else {
      gmate->AddElement(el, wmat[i]);
    }
  }
  // add the material to the List
  G3Mat.put(&imate,gmate);
}


