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
/// \file parameterisations/gflash/src/ExGflashMaterialManager.cc
/// \brief Implementation of the ExGflashMaterialManager class
//
#include "ExGflashMaterialManager.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

ExGflashMaterialManager* ExGflashMaterialManager::mpointer=0;

void ExGflashMaterialManager::storeElement(G4String name,G4String symbol, double z, double a)
{
  if (elist.find(symbol)!=elist.end())
  {
    G4cout<<"attempt of redefining an existing element "<<symbol<<G4endl;
                G4cout<<elist[symbol]<<G4endl;      
  }
  else
  {
     G4Element *newelement=new G4Element(name,symbol,z,a);
    elist[symbol]=newelement;
  }
}

G4Element* ExGflashMaterialManager::getElement(G4String symbol)
{
  if (elist.find(symbol)==elist.end())
  {
    G4cout<<" element "<<symbol<<" not found!"<<G4endl;
    return NULL;
  }
  else
    return elist[symbol];
}

void ExGflashMaterialManager::storeMaterial(G4String name, double z, double a, double d)
{
  if (mlist.find(name)!=mlist.end())
  {
    G4cout<<"attempt of redefining an existing material "<<G4endl;
                G4cout<<mlist[name]<<G4endl;      
  }
  else
  {
     G4Material *newmaterial=new G4Material(name,z,a,d);
    mlist[name]=newmaterial;
  }
}

void ExGflashMaterialManager::storeMaterial(G4String name, double d, int ncomponent)
{
  if (mlist.find(name)!=mlist.end())
  {
    G4cout<<"attempt of redefining an existing material "<<G4endl;
                G4cout<<mlist[name]<<G4endl;      
  }
  else
  {
     G4Material *newmaterial=new G4Material(name,d,ncomponent);
    mlist[name]=newmaterial;
  }
}

G4Material* ExGflashMaterialManager::getMaterial(G4String name)
{
  if (mlist.find(name)==mlist.end())
  {
    G4cout<<" material "<<name<<" not found!"<<G4endl;
    return NULL;
  }
  else
  {
    // G4cout<<"returning material "<<name<<G4endl;
    return mlist[name];
  }
}

void ExGflashMaterialManager::addMaterial(G4String n1,G4String n2,double fraction)
{
  G4Material *mat=getMaterial(n1);
  mat->AddMaterial(getMaterial(n2),fraction);
}

void ExGflashMaterialManager::addElement(G4String n1,G4String n2,double fraction)
{
  G4Material *mat=getMaterial(n1);
  mat->AddElement(getElement(n2),fraction);
}

void ExGflashMaterialManager::addElement(G4String n1,G4String n2,int natoms)
{
  G4Material *mat=getMaterial(n1);
  mat->AddElement(getElement(n2),natoms);
}

void ExGflashMaterialManager::printMaterialTable()
{
  MaterialList::iterator it;
  for (it=mlist.begin();it!=mlist.end();it++) G4cout<<(*it).second<<G4endl;  
}

void ExGflashMaterialManager::printElementTable()
{
  ElementList::iterator it;
  for (it=elist.begin();it!=elist.end();it++) G4cout<<(*it).second<<G4endl;
}

void ExGflashMaterialManager::initialize()
{
        G4String name,symbol;
        double z,a,density;
        int ncomponents,natoms;
  //double fraction;
  storeElement(name="Hydrogen",symbol="H" , z= 1., a=1.01*g/mole);
  storeElement(name="Carbon"  ,symbol="C" , z= 6., a=12.01*g/mole);
  storeElement(name="Nitrogen",symbol="N" , z= 7., a=14.01*g/mole);
  storeElement(name="Oxygen"  ,symbol="O" , z= 8., a=16.00*g/mole);
  storeElement(name="Silicon",symbol="Si" , z= 14., a=28.09*g/mole);
  storeElement(name="Argon",symbol="Ar",z=18.,a=39.95*g/mole);
  storeElement(name="Iron"    ,symbol="Fe", z=26., a=55.85*g/mole);
  storeElement(name="Aluminum",symbol="Al",z=13.,a=26.98*g/mole);
  storeElement(name="Lead",symbol="Pb",z=82.,a=207.19*g/mole);
  storeElement(name="Fluorine",symbol="F",z=9.,a=18.99*g/mole);
  storeElement(name="Chlorine",symbol="Cl",z=17.,a=35.45*g/mole);
        storeElement(name="Tungsten",symbol="W",z=74.,a=183.85*g/mole);
  storeMaterial(name="Aluminium", z=13., a = 26.98*g/mole, density = 2.700*g/cm3);
  storeMaterial(name="Iron",z=26., a=55.85*g/mole, density=7.87*g/cm3);
  storeMaterial(name="Copper",z=29.,a=63.546*g/mole, density=8.96*g/cm3);
  storeMaterial(name="Silicon",z=14.,a=28.0855*g/mole, density=2.33*g/cm3);
  storeMaterial(name="Tungsten",z=74.,a=183.85*g/mole,  density = 19.3*g/cm3 );
  storeMaterial(name="Lead"     , z=82., a= 207.19*g/mole,density = 11.35*g/cm3);
  storeMaterial(name="LAr",z=18.,a=39.95*g/mole,density =1.39*g/cm3);
  storeMaterial(name="Scintillator", density=1.032*g/cm3, ncomponents=2);
  addElement("Scintillator","C",natoms=9);
  addElement("Scintillator","H",natoms=10);
  storeMaterial(name="PbWO4", density=8.28*g/cm3, ncomponents=3);
   addElement("PbWO4","Pb",natoms=1);
   addElement("PbWO4","W",natoms=1);
   addElement("PbWO4","O",natoms=4);
  //  addElement("PbWO4","Pb", fraction=0.45532661 );
  //  addElement("PbWO4","W",  fraction=0.40403397 );
  //  addElement("PbWO4","O",  fraction=0.14063942);
  storeMaterial(name="Air"  , density=1.290*mg/cm3, ncomponents=2);
  addElement("Air", "N", .7);
  addElement("Air", "O", .3);
  storeMaterial(name="CO2", density=1.977*mg/1000*cm3,ncomponents=2);
  addElement("CO2","C",natoms=1);
  addElement("CO2","O",natoms=2);
  storeMaterial(name="ArCO2",density=1.8*mg/1000*cm3,ncomponents=2);
  addElement("ArCO2","Ar",.93);
  addMaterial("ArCO2","CO2",.07);
  storeMaterial(name="RPCgas", density=1.977*mg/1000*cm3,ncomponents=3);
  addElement("RPCgas","C",natoms=2);
  addElement("RPCgas","H",natoms=2);
  addElement("RPCgas","F",natoms=4);
  storeMaterial(name="RPVC", density=1.4*g/cm3,ncomponents=3);
  addElement("RPVC","C",natoms=2);
  addElement("RPVC","H",natoms=3);
  addElement("RPVC","Cl",natoms=1);
  storeMaterial(name="Bakelite", density=1.4*g/cm3,ncomponents=3);
  addElement("Bakelite","C",natoms=1);
  addElement("Bakelite","H",natoms=4);
  addElement("Bakelite","O",natoms=2);
}












