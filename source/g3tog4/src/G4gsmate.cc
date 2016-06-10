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
// $Id: G4gsmate.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, 27 Sep 99

#include <cmath>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3EleTable.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"

void PG4gsmate(G4String *tokens)
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
  // G4double absl = Rpar[4];
  G4double *ubuf = &Rpar[5];

  G4gsmate(imate, name, a, z, dens, radl, nwbf, ubuf);
}

/*
// replaced with G3EleTable 
// only used G4Elements are created;
// !! no checking of given A of the element;
//

G4Element* CreateElement(G4double zeff, G4double aeff, G4String matName)
{
  // tolerance in Z, A for element definition
  const G4double tolerance = 0.001;

  // define the symbol Z%A% of element
  // !! symbol is not unambiguous element identifier
  char symbol[20];
  sprintf(symbol,"Z%dA%d",int(zeff),int(aeff/(g/mole)));
  G4String elSymbol = symbol;
 
  // search element table for the element with given (Z,A)
  // 
  G4int index = 0;
  const G4ElementTable* table = G4Element::GetElementTable();
  for (G4int i=0; i<table->entries(); i++) {
     G4Element* entry = (*table)[i];
     if (elSymbol == entry->GetSymbol()) index++; 
     if ( std::abs(zeff - entry->GetZ()) < tolerance &&
         (std::abs(aeff - entry->GetA())/(g/mole)) < tolerance ){
       return entry;
     }
  }  

  // define a unique name En-Z%A% 
  // (n - index of elements with the same int(Z), int(A))
  char chIndex[4];
  sprintf(chIndex,"%d",index);
  G4String elName = "E";  
  elName = elName + chIndex + "-";
  elName = elName + elSymbol;  

  // create new element if it was not found in element table
  G4Element* element = new G4Element(elName, elSymbol, zeff, aeff);	
  G4cout << "New element: " << element->GetName()
         << " for " << matName << " material has been created." << G4endl;
  return element;	
}
*/

void G4gsmate(G4int imate, G4String name, G4double ain, G4double zin,
              G4double densin, G4double, G4int, G4double*)
{
  G4double G3_minimum_density = 1.e-10*g/cm3;

  // add units
  G4double z = zin;    
  G4double a = ain*g/mole;
  G4double dens = densin*g/cm3;

  G4Material* material=0;
  
  G4String sname = name.strip(G4String::both);
  if (sname == "AIR") {
    // handle the built in AIR mixture
    G4double aa[2], zz[2], wmat[2];
    aa[0] = 14.01*g/mole;
    aa[1] = 16.00*g/mole;
    zz[0] = 7;
    zz[1] = 8;
    wmat[0] = 0.7;
    wmat[1] = 0.3;
    // G4double theDensity = 1.2931*mg/cm3;
    G4double theDensity = 0.0012931;
    G4int n=2;
    G4gsmixt(imate, sname, aa, zz, theDensity, n, wmat);
  } 
  else if ( z<1 || dens < G3_minimum_density ) {
    // define vacuum according to definition from N03 example
    G4double density     = universe_mean_density;    //from PhysicalConstants.h
    G4double pressure    = 3.e-18*pascal;
    G4double temperature = 2.73*kelvin;
    material = new G4Material(name, z=1., a=1.01*g/mole, density,
                    kStateGas,temperature,pressure);
  }
  else {
    //G4Element* element = CreateElement(z, a, name);
    G4Element* element = G3Ele.GetEle(z);
    material = new G4Material(name, dens, 1);
    material->AddElement(element, 1.);    
  }  

  // add the material to the List
  G3Mat.put(imate, material);
}







