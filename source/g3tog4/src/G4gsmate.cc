// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsmate.cc,v 1.4 1999-12-05 17:50:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include <math.h>

#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3EleTable.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"

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
     if ( abs(zeff - entry->GetZ()) < tolerance &&
         (abs(aeff - entry->GetA())/(g/mole)) < tolerance ){
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
         << " for " << matName << " material has been created." << endl;
  return element;	
}
*/

void G4gsmate(G4int imate, G4String name, G4double ain, G4double zin,
              G4double densin, G4double radl, G4int nwbf, G4double* ubuf)
{
  G4double G3_minimum_density = 1.e-10*g/cm3;

  // add units
  G4double z = zin;    
  G4double a = ain*g/mole;
  G4double dens = densin*g/cm3;

  G4Material* material;
  
  G4String sname = name.strip(G4String::both);
  if (sname == "AIR") {
    // handle the built in AIR mixture
    G4double a[2], z[2], wmat[2];
    a[0] = 14.01*g/mole;
    a[1] = 16.00*g/mole;
    z[0] = 7;
    z[1] = 8;
    wmat[0] = 0.7;
    wmat[1] = 0.3;
    // G4double theDensity = 1.2931*mg/cm3;
    G4double theDensity = 0.0012931;
    int n=2;
    G4gsmixt(imate, sname, a, z, theDensity, n, wmat);
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







