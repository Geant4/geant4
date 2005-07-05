//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst51DetectorConstruction.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// author: Susanna Guatelli (guatelli@ge.infn.it)
// 
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#include "Tst51DetectorConstruction.hh"
#include "Tst51DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4TransportationManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"

Tst51DetectorConstruction::Tst51DetectorConstruction()
  :hydrogen(0),beryllium(0), 
   aluminium(0),silicon(0),iron(0), germanium(0),
   gold(0), lead(0), air(0),vacuum(0),
   targetMaterial(0),defaultMaterial(0),
   solidWorld(0),logicWorld(0),physiWorld(0),
   solidTarget(0),logicTarget(0),physiTarget(0)
{
  // default parameter values of the target
  targetThickness = 0.0094 * cm;
  targetY = 20.*m;
  targetX = 20.*m;
 
  messenger = new Tst51DetectorMessenger(this);
}

Tst51DetectorConstruction::~Tst51DetectorConstruction()
{
  delete messenger;
  delete vacuum;
  delete air;
  delete lead;
  delete gold;
  delete germanium;
  delete iron; 
  delete silicon;
  delete aluminium;
  delete beryllium;
  delete hydrogen;
}

G4VPhysicalVolume* Tst51DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructWorld();
}

void Tst51DetectorConstruction::DefineMaterials()
{ 
 
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  
 
  G4int ncomponents;
  G4double fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  // Materials

  a = 1.01*g/mole;
  density = 8.3748e-5 *g/cm3; 
  hydrogen = new G4Material(name="Hydrogen", z=1., a, density);

 
  a = 9.012*g/mole;
  density = 1.848*g/cm3; 
  beryllium = new G4Material(name="Beryllium", z=4., a, density);


  a = 16.00*g/mole;
  density = 0.00133*g/cm3;
  ossigeno = new G4Material(name="Oxygen", z=8.,a,density);


  a = 26.981*g/mole;
  density = 2.6989*g/cm3; 
  aluminium = new G4Material(name="Aluminium", z=13., a, density);
 

  a = 28.085*g/mole;
  density = 2.33*g/cm3; 
  silicon = new G4Material(name="Silicon", z=14., a, density);


  a = 55.845*g/mole;
  density = 7.874*g/cm3;
  iron = new G4Material(name="Iron", z=26., a, density);

 
  density = 5.323*g/cm3;
  a = 72.64*g/mole;
  germanium = new G4Material(name="Germanium", z=32., a, density);


  density = 19.32*g/cm3;
  a = 196.966*g/mole;
  gold = new G4Material(name="Gold", z=79., a, density);

 
  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  lead = new G4Material(name="Lead", z=82., a, density);

 
  density = 1.290*mg/cm3;
  air = new G4Material(name="Air"  , density, ncomponents=2);
  air->AddElement(N, fractionmass=0.7);
  air->AddElement(O, fractionmass=0.3);

  //
  // examples of vacuum
  //

  density = universe_mean_density;
  pressure = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole,
				      density,kStateGas,temperature,pressure);

  targetMaterial = aluminium;
  defaultMaterial  = vacuum;
}
  
G4VPhysicalVolume* Tst51DetectorConstruction::ConstructWorld()
{
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
			 50.*m,50.*m,50.*m);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
 
  // Slab filled with absorber material ...
  solidTarget = 0; 
  logicTarget = 0; 
  physiTarget = 0;  
  
  if (targetThickness > 0.) 
    {solidTarget = new G4Box("Target",		//its name
			      targetX/2,targetY/2, targetThickness/2 ); 
                          
    logicTarget = new G4LogicalVolume(solidTarget,    //its solid
				      targetMaterial, //its material
				      "Target");      //its name
      			                  
    }

    physiTarget = new G4PVPlacement(0,		   //no rotation
				    G4ThreeVector(0.,0.,0.),  //its position
				    "Target",        //its name
				    logicTarget,     //its logical volume
				    physiWorld,        //its mother
				    false,             //no boulean operat
				    0);                //copy number 

  G4Colour  magenta (1.0, 0.0, 1.0) ; 

  G4VisAttributes* simpleVisAtt= new G4VisAttributes(magenta);
  simpleVisAtt->SetVisibility(true);
  simpleVisAtt->SetForceSolid(true);
  logicTarget->SetVisAttributes(simpleVisAtt);  

  PrintParameters();
  return physiWorld;
    }

void Tst51DetectorConstruction::PrintParameters()
{
  G4cout << " Target z size is: "
        
         << targetThickness/cm<<" cm" <<G4endl; 

  G4cout<<" Target y size is: "<<targetY/cm<<" cm"<<G4endl;
  G4cout<<" Target z size is: "<<targetX/cm<<" cm"<<G4endl;
  G4cout<<" Target material is: "<<targetMaterial->GetName()<<G4endl; 
}

void Tst51DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
    {
    targetMaterial = pttoMaterial;
    logicTarget->SetMaterial(pttoMaterial); 
    PrintParameters();
    }
  else {G4cout<<materialChoice <<"is not available!!!!"<< G4endl;}  
  UpdateGeometry();           
}

void Tst51DetectorConstruction::SetTargetThickness(G4double val)
{ 
 
  // change Target thickness
  targetThickness = val;
  UpdateGeometry();
  G4cout<<"The thickness of the target is: "<< targetThickness/cm << " cm"<<G4endl;
}  
G4double Tst51DetectorConstruction::GetTargetThickness()
{ 
  // Get target thickness
  return  targetThickness;
}  

void Tst51DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructWorld());
}

