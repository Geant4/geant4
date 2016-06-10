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
// $Id: G4gsmixt.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, 27 Sep 99

#include <iomanip>
#include <iomanip>

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G3toG4.hh"
#include "G3EleTable.hh"
#include "G3MatTable.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"

void PG4gsmixt(G4String *tokens)
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
    G4double *z = Rpar + 1+std::abs(nlmat);
    G4double *wmat = Rpar + 1 + 2*std::abs(nlmat);
/*
    for (int i=0; i<std::abs(nlmat); i++){
      //Rpar[i]=Rpar[i]*g/mole;
      Rpar[i]=Rpar[i];
    };
*/
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
    for (i=0; i<std::abs(nlmat); i++) { 
      // total molecular weight 
      aMol += wmat[i]*a[i]; 
    }  
    if (aMol == 0.) {
      G4String text = "G4mixt: Total molecular weight in " + name + " = 0.";       
      G4Exception("G4gsmixt()", "G3toG40016", FatalException, text);
      return;
    }
    for (i=0; i<std::abs(nlmat); i++) {
      // weight fractions
      wmat[i] = wmat[i]*a[i]/aMol;
    }
  }

  // create material with given number of components
  // (elements)

  G4Material* material 
    = new G4Material(name, dens*g/cm3, std::abs(nlmat));
  for (i=0; i<std::abs(nlmat); i++) {
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
  G4int nmate = std::abs(nlmat);
  G4String sname = name.strip(G4String::both);
  G4double theDensity = dens*g/cm3;

  G4Material* theMixture = new G4Material(name, dens, nmate); 
  G4bool ok=true;
  for (int i=0; i< nmate; i++){
    G4Element* theElement = G3Ele.GetEle(z[i]);
    if (nlmat>0) {
      G4double fractionmass = wmat[i];
      ok = ok && std::abs(fractionmass)<=1.;
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
	     << "' some |weights|>1 : " << G4endl;
      for (G4int i=0;i<nlmat; i++) {
	G4cerr << "Component " << std::setw(3) << i+1 << " fraction: "
	       << std::setw(10) << wmat[i] << G4endl;
      }
    } else if (nlmat<0) {
      G4cerr << "G4gsmixt: for mixture '" << name 
	     << "' some #natoms are non-integer: " << G4endl;
      for (G4int i=0;i<nlmat; i++) {
	G4cerr << "Component " << std::setw(3) << i+1 << " #atoms "
	       << std::setw(10) << wmat[i] << G4endl;
      }
    } else {
      G4cerr << "G4gsmixt: Number of components for mixture '" 
	     << name << "' (" << nlmat << ") not allowed." << G4endl;
    }
  }
}
*/



