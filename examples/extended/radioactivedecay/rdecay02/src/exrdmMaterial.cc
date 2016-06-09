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
/// \file radioactivedecay/rdecay02/src/exrdmMaterial.cc
/// \brief Implementation of the exrdmMaterial class
//
////////////////////////////////////////////////////////////////////////////////
//
#include "exrdmMaterial.hh"
#include "exrdmMaterialData.hh"
#include "exrdmMaterialMessenger.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <vector>
#include <iomanip>  
////////////////////////////////////////////////////////////////////////////////
//
exrdmMaterial::exrdmMaterial ()
{
  fMaterial.clear();
  fElement.clear();
  fIsotope.clear();
  // some default materials vacuum (0), air (1)  and aluminium (2)  defined here
  // examples of vacuum
  //
  //  G4double a,z;

  static G4bool bmat = false ;

  if (!bmat) {

    // vacuum
    G4double  density    = universe_mean_density;    //from PhysicalConstants.h
    G4double pressure    = 3.e-18*pascal;
    G4double temperature = 2.73*kelvin;
    AddMaterial("Vacuum", "H", density,"gas",temperature,pressure);
    // air
    density   = 1.290*mg/cm3;
    AddMaterial("Air", "N0.78-O0.22", density, "gas");

    // aluminium
    density=2.700*g/cm3 ;
    AddMaterial ("Aluminium", "Al", density,"");
    
     //silicon
    density=2.3290*g/cm3 ;
    AddMaterial ("Silicon", "Si", density,"");

    bmat              = true;
  }
  // create commands for interactive definition of the geometry 
  fMaterialMessenger = new exrdmMaterialMessenger(this);
}
////////////////////////////////////////////////////////////////////////////////
//
exrdmMaterial::~exrdmMaterial ()
{
  delete fMaterialMessenger;
}
////////////////////////////////////////////////////////////////////////////////
//
void exrdmMaterial::AddMaterial (G4String name, G4String formula,
                G4double density, G4String state, G4double tem, G4double pres)
{
  G4int isotope, Z;
  size_t i;
  for (i = 0; i<fMaterial.size(); i++) {
    if (fMaterial[i]->GetName() == name) {
      G4cerr <<" AddMaterial : material " <<name
             <<" already exists." <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
      return;
    }
  }
  
  char *tokenPtr1 = NULL;
  char *sname     = NULL;
  G4String s0, s1("0123456789");
  G4String element, isotopename;
  G4int ncomponents, natoms;
  G4double fatoms = 0.;
  size_t ls, id=0, ll, lr;
  ncomponents = 0;

  sname       = new char[strlen(formula)+1];
  strcpy(sname,formula);
  tokenPtr1 = strtok(sname,"-");

  while (tokenPtr1 != NULL) {
    ncomponents++;
    tokenPtr1 = strtok( NULL, "-");
  }
  delete[] sname;

  G4Material* aMaterial = 0;
  G4cout << name <<" "<< formula << " " << density/(g/cm3) << " " << tem <<" "
                                                                       <<pres << G4endl;
 
  if (state == "") {
    aMaterial = new G4Material(name, density, ncomponents);
  } else if (state == "solid" && tem > 0.) {
    aMaterial = new G4Material(name, density, ncomponents, 
                                           kStateSolid, tem );
  } else if (state == "gas" && pres > 0.) {
    aMaterial = new G4Material(name, density, ncomponents, 
                                           kStateGas, tem, pres );
  }
  if (aMaterial == 0) {
    G4cerr <<" AddMaterial : Name " <<name <<"." <<G4endl;
    G4cerr <<"--> Command failed." <<G4endl;
    return;
  }

  sname=new char[strlen(formula)+1];
  strcpy(sname,formula);
  tokenPtr1 = strtok(sname,"-");

  while (tokenPtr1 != NULL) {
    isotope = 0;
    //      G4cout << tokenPtr1 << G4endl;
    s0      = G4String(tokenPtr1);
    ls      = s0.length();
    ll      = s0.find("(");
    lr      = s0.find(")");
    if (ll == lr) {
      id = s0.find_first_of(s1);
      element = s0.substr(0,id);
      
      if (element.length() == 1) element.insert(0," ");
      for (i = 0; i<110; i++) {
        if (element == fELU[i]) break;
      }
      if (i == 110) {
        for (i = 0; i<110; i++) {
          if (element == fELL[i]) break;
        }
        if (i == 110) {
          for (i = 0; i<110; i++) {
            if (element == fEUU[i]) break;
          }
        }
      }
      
      if (i == 110) {
        G4cerr <<"AddMaterial : Invalid element in material formula."
               <<element <<G4endl;
        G4cerr <<"--> Command rejected." <<G4endl;
//        delete aMaterial;
//        fMaterial[NbMat] = NULL;
        return;
      }

      Z       = i+1;
      element = fELU[i];
      if (id == std::string::npos) {
        natoms = 1;
      } else {
        natoms = atoi((s0.substr(id,ls-id)).c_str());
      }
      if (natoms < 1) fatoms = atof((s0.substr(id,ls-id)).c_str());
      //        G4cout << "   Elements = " << element << G4endl;
      //G4cout << "   Nb of atoms = " << natoms << G4endl;
    } else {
      element = s0.substr(0,ll);
      isotope = atoi((s0.substr(ll+1,lr-ll)).c_str());
      if (element.length() == 1) element.insert(0," ");
      for (i = 0; i<110; i++) {
        if (element == fELU[i]) break;
      }
      if (i == 110) {
        for (i = 0; i<110; i++) {
          if (element == fELL[i]) break;
        }
        if (i == 110) {
          for (i = 0; i<110; i++) {
            if (element == fEUU[i]) break;
          }
        }
      }
      if (i == 110) {
        G4cerr <<"AddMaterial : Invalid element in material formula."
               <<element <<G4endl;
        G4cerr <<"--> Command rejected." <<G4endl;
//        delete aMaterial;
//        fMaterial[NbMat] = NULL;
        return;
      }

      Z           = i+1;
      element     = fELU[i];
      isotopename = element+s0.substr(ll+1,lr-ll-1);
      if (lr == std::string::npos ) {
        natoms = 1;
      } else {
        natoms = atoi((s0.substr(lr+1,ls-lr)).c_str());
      }  
      if (natoms < 1)  fatoms = atof((s0.substr(id,ls-id)).c_str());
      if (fatoms == 0.) natoms = 1;
      //
      //        G4cout << "   Elements = " << element << G4endl;
      //   G4cout << "   fIsotope Nb = " << isotope << G4endl;
      //        G4cout << "   Nb of atoms = " << natoms << G4endl;
    }
    if (isotope != 0) {
      if (G4Isotope::GetIsotope(isotopename) == NULL) {
        G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope, isotope*g/mole);
        G4Element* aElement = new G4Element(isotopename, element, 1);
        aElement->AddIsotope(aIsotope, 100.*perCent);
        fIsotope.push_back(aIsotope);
        if (natoms>0) { 
          aMaterial->AddElement(aElement, natoms);
        } else {
          aMaterial->AddElement(aElement, fatoms);
        }
        fElement.push_back(aElement);
      } else {
        if (natoms>0) { 
          aMaterial->AddElement( G4Element::GetElement(isotopename,false) , natoms);
        } else {
          aMaterial->AddElement( G4Element::GetElement(isotopename,false) , fatoms);
        }
      }      
    } else {
      if ( G4Element::GetElement(element,false) == NULL) {
        G4Element* aElement = new G4Element(element, element, Z, fA[Z-1]*g/mole);
        if (natoms>0) { 
          aMaterial->AddElement(aElement, natoms);
        } else {
          aMaterial->AddElement(aElement, fatoms);
        }
        fElement.push_back(aElement);
      } else {
        if (natoms>0) { 
          aMaterial->AddElement( G4Element::GetElement(element,false) , natoms);
        } else {
          aMaterial->AddElement( G4Element::GetElement(element,false) , fatoms);
        } 
      }
    }
    tokenPtr1 = strtok( NULL, "-");
    //      s0.empty();
    //element.erase();
    //
  }

  delete[] sname;
  fMaterial.push_back(aMaterial);
  G4cout <<" fMaterial:" <<name <<" with formula: " <<formula <<" added! "
         <<G4endl;
  G4cout <<"     Nb of fMaterial = " <<fMaterial.size() <<G4endl;
  G4cout <<"     Nb of fIsotope =  " <<fIsotope.size() <<G4endl;
  G4cout <<"     Nb of fElement =  " <<fElement.size() <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void exrdmMaterial::DeleteMaterial (G4int j)
{
  size_t i(j-1);
  if (i > fMaterial.size()) {
    G4cerr <<"DeleteMaterial : Invalid material index " <<j <<"." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void exrdmMaterial::DeleteMaterial (G4String )
{
  G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
         <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int exrdmMaterial::GetMaterialIndex (G4String name)
{
  size_t i ;
  for (i = 0; i < fMaterial.size(); i++) {
    if (fMaterial[i]->GetName() == name) break;
  }
  G4int k = G4int(i);
  if (i == fMaterial.size()) k = -1;
  return k;
}
////////////////////////////////////////////////////////////////////////////////
//
void exrdmMaterial::ListMaterial ()
{
  G4cout <<" There are" <<std::setw(3) <<fMaterial.size()
         <<" materials defined."  <<G4endl;
  for (size_t i = 0; i< fMaterial.size(); i++) 
    G4cout <<"     fMaterial Index " <<std::setw(3) <<i+1 <<" "
           <<std::setw(14) <<fMaterial[i]->GetName()
           <<"  density: " <<std::setw(6) <<std::setprecision(3)
           <<G4BestUnit(fMaterial[i]->GetDensity(),"Volumic Mass") <<G4endl;
  G4cout <<G4endl;
  
}
////////////////////////////////////////////////////////////////////////////////
