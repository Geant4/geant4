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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include <fstream>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyModulator.hh"
#include "G4Transform3D.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include <iostream>

using namespace std;

HadrontherapyModulator::HadrontherapyModulator(): physiMotherMod(0),
						 solidMod1(0),         logicMod1(0),          physiMod1(0),
						 solidMod2(0),         logicMod2(0),          physiMod2(0),
						 solidMod3(0),         logicMod3(0),          physiMod3(0),
						 solidMod4(0),         logicMod4(0),          physiMod4(0),
						 FileName("Modulators/Modulator009_02.txt")
{ 
   pi=4*atan(1.);
   StepNumbers=22;
   Weight=new G4double[StepNumbers];
   StepThickness=new G4double[StepNumbers];
   StartingAngle=new G4double[StepNumbers];
   SpanningAngle=new G4double[StepNumbers];
   PositionMod=new G4ThreeVector[StepNumbers];


   solidMod=new G4Tubs *[StepNumbers];
   logicMod=new G4LogicalVolume *[StepNumbers];
   physiMod=new G4VPhysicalVolume *[(4*(StepNumbers-1)+1)];
     
   for (G4int i=0;i<StepNumbers;i++)
  {
	Weight[i]=0;
	StepThickness[i]=0;
	StartingAngle[i]=0;
	SpanningAngle[i]=0;
	PositionMod[i]=G4ThreeVector(0,0,0);
	solidMod[i]=0;
	logicMod[i]=0;  
	  
  }
  
  for (G4int i=0;i<4*(StepNumbers-1)+1;i++)
  {
  physiMod[i]=0;	  
  }
	
  
  ModulatorMessenger = new  HadrontherapyModulatorMessenger(this);	
  ModulatorDefaultProperties();
  rm = new G4RotationMatrix(); 
  G4double phi = 270. *deg;     
  rm -> rotateY(phi); 
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyModulator::~HadrontherapyModulator() 
{
  delete rm;
  delete [] Weight;
  delete []	StepThickness;
  delete []	StartingAngle;
  delete []	SpanningAngle;
  delete []	PositionMod;
  delete []	solidMod;
  delete []	logicMod;  
  delete []    physiMod;
  delete ModulatorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator::ModulatorDefaultProperties()
{
/* Here we initialize the step properties of Modulator wheel, you can create your
specific modulator by changing the values in this class or writing them in an external
file and activate reading from file via a macrofile.	*/	

 StepThickness[0]=0; Weight[0]=.14445;	
 StepThickness[1]=.8; Weight[1]=.05665;	
 StepThickness[2]=1.6; Weight[2]=.05049;
 StepThickness[3]=2.4; Weight[3]=.04239;
 StepThickness[4]=3.2; Weight[4]=.04313;
 StepThickness[5]=4.0; Weight[5]=.03879;
 StepThickness[6]=4.8; Weight[6]=.04182;
 StepThickness[7]=5.6; Weight[7]=.03422;
 StepThickness[8]=6.4; Weight[8]=.03469;
 StepThickness[9]=7.2; Weight[9]=.03589;
 StepThickness[10]=8.0; Weight[10]=.03633;
 StepThickness[11]=8.8; Weight[11]=.03842;
 StepThickness[12]=9.6; Weight[12]=.03688;
 StepThickness[13]=10.4; Weight[13]=.03705;
 StepThickness[14]=11.2; Weight[14]=.03773;
 StepThickness[15]=12.0; Weight[15]=.03968;
 StepThickness[16]=12.8; Weight[16]=.04058;
 StepThickness[17]=13.6; Weight[17]=.03903;
 StepThickness[18]=14.4; Weight[18]=.04370;
 StepThickness[19]=15.2; Weight[19]=.03981;
 StepThickness[20]=16.0; Weight[20]=.05226;
 StepThickness[21]=16.8; Weight[21]=.03603;
 GetStepInformation();	
 
} 
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator:: ModulatorPropertiesFromFile(G4String Name)
{
  delete [] Weight;
  delete []	StepThickness;
  delete []	StartingAngle;
  delete []	SpanningAngle;
  delete []	PositionMod;
  delete []	solidMod;
  delete []	logicMod;  
  delete []    physiMod;
  delete    solidMod1;
  delete    logicMod1;
  delete    physiMod1;
  delete    solidMod2;       
  delete    logicMod2;          
  delete    physiMod2;
  delete    solidMod3;         
  delete    logicMod3;       
  delete    physiMod3;
  delete    solidMod4;       
  delete    logicMod4;       
  delete    physiMod4;
// The Modulator wheel properties is getting form an external file "ModoulWeight.txt" 
  File.open(Name,  std::ios::in);
  if(!File.is_open())
  {
  G4cout<<" WARNING: The File with name of "<<Name<<
 " doesn't exist to get modulator step properties. please modify it and try again"<<G4endl;
 
 G4Exception("HadrontherapyModulator::ModulatorPropertiesFromFile( )", "Hadrontherapy0009"
 , FatalException, "Error: No available external file for reading from");
    }	  

  G4String string;
  File >>string>> StepNumbers;
  File >>string>>string>>string;
  
  
   Weight=new G4double[StepNumbers];
   StepThickness=new G4double[StepNumbers];
   StartingAngle=new G4double[StepNumbers];
   SpanningAngle=new G4double[StepNumbers];
   PositionMod=new G4ThreeVector[StepNumbers];


   solidMod=new G4Tubs *[StepNumbers];
   logicMod=new G4LogicalVolume *[StepNumbers];
   physiMod=new G4VPhysicalVolume *[(4*(StepNumbers-1)+1)];
 
  for(G4int i=0;i<StepNumbers;i++)
   {
	 G4String stringX;
	 File>>stringX>> StepThickness[i]>>Weight[i];  
   }	
   
   GetStepInformation();
   BuildSteps();
   


}
////////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator::GetStepInformation()
{

 G4double TotalWeight=0;
 // convert the absolute weight values to relative ones 
 for(G4int i=0;i<StepNumbers;i++)
  {
	 TotalWeight+=Weight[i];
  }
 
 for(G4int i=0;i<StepNumbers;i++)
 {
   Weight[i]=Weight[i]/TotalWeight;
 }
 
 // To build the RMW step layers will be put one after another  
 
  StartingAngle[0]=0 *deg;
  SpanningAngle[0]=90 *deg;
  G4double PositionModx;
  G4double WholeStartingAngle=0 *deg;
  G4double WholeThickness=0;
  for(G4int i=1;i<StepNumbers;i++)
  {
	  StartingAngle[i]=WholeStartingAngle+(Weight[i-1]*(2*pi))/8;
	  SpanningAngle[i]=90* deg -2*StartingAngle[i];
	  StepThickness[i]=StepThickness[i]-WholeThickness;
	  PositionModx=WholeThickness+StepThickness[i]/2.;
	  PositionMod[i]=G4ThreeVector(0,0,PositionModx);
	  WholeThickness+=StepThickness[i];
	  WholeStartingAngle=StartingAngle[i];
  }	
	
	
}
/////////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator::BuildModulator(G4VPhysicalVolume* motherVolume)
{
  G4bool isotopes = false;
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);


   Mod0Mater = airNist;
   ModMater = airNist; // You have to change modulator material via a macrofile (default is air)
 
  innerRadiusOfTheTube = 2.5 *cm;
  outerRadiusOfTheTube = 9.5 *cm;

  // Mother of the modulator wheel  
  G4ThreeVector positionMotherMod = G4ThreeVector(-2160.50 *mm, 30 *mm, 50 *mm);
 
  G4Box* solidMotherMod = new G4Box("MotherMod", 12 *cm, 12 *cm, 12 *cm);
 
  logicMotherMod = new G4LogicalVolume(solidMotherMod, Mod0Mater,"MotherMod",0,0,0);

  physiMotherMod = new G4PVPlacement(rm,positionMotherMod,  "MotherMod", 
				     logicMotherMod,    				  
				     motherVolume,      
				     false,           
				     0); 
	BuildSteps();			      
				     
				     
				     
				 }            
 ///////////////////////////////////////////////////////////////////////////////////////
 void HadrontherapyModulator::BuildSteps()
 {
  //----------------------------------------------------------
  // Mother volume of first quarter of the modulator
  //----------------------------------------------------------

  G4double hightOfTheTube0 = 10.0 *cm;
  G4double startAngleOfTheTube0 = 0 *deg;
  G4double spanningAngleOfTheTube0 = 90 *deg;
  
  G4RotationMatrix rm1;
  rm1.rotateZ(0 *deg);
 
  G4ThreeVector positionMod1 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod1 = new G4Tubs("Mod1",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube0/2.,
			 startAngleOfTheTube0, 
			 spanningAngleOfTheTube0);
  
  logicMod1 = new G4LogicalVolume(solidMod1, Mod0Mater, "Mod1",0,0,0);
  
  physiMod1 = new G4PVPlacement(G4Transform3D(rm1, positionMod1), 
				logicMod1,    
				"Mod1",       
				logicMotherMod,  
				false,         
				0);            
  
 
  //----------------------------------------------------------
  //  modulator steps
  //----------------------------------------------------------
  for (G4int i=1;i<StepNumbers;i++)
  {

  solidMod[i] = new G4Tubs("Modstep",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 StepThickness[i]/2.,
			 StartingAngle[i], 
			 SpanningAngle[i]);
			       
  logicMod[i] = new G4LogicalVolume(solidMod[i], 
                                   ModMater, "Modstep",0,0,0);
  
  physiMod[i] = new G4PVPlacement(0,               
				PositionMod[i],  
				logicMod[i],     
				"Modstep",        
				logicMod1,     
				false,         
				0);     
	  
	  
  }
 
 /////////////////////////////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------
  // Mother volume of the second modulator quarter
  //----------------------------------------------------------
  
    
  G4RotationMatrix rm2;
  rm2.rotateZ(90 *deg);
  
  G4ThreeVector positionMod2 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod2 = new G4Tubs("Mod2",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0/2.,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  logicMod2 = new G4LogicalVolume(solidMod2,
                                  Mod0Mater, "Mod2",0,0,0);
  
  
  physiMod2 = new G4PVPlacement(G4Transform3D(rm2, positionMod2), 
				 logicMod2,    
				 "Mod2",       
				 logicMotherMod, 
				 false,         
				 0);            
 

    for (G4int i=1;i<StepNumbers;i++)
  {
	
  physiMod[StepNumbers+i-1] = new G4PVPlacement(0,               
				PositionMod[i],  
				logicMod[i],     
				"Modstep",        
				logicMod2,     
				false,         
				0);  
				 
			} 

 

  //----------------------------------------------------------
  // Mother volume of the third modulator quarter
  //----------------------------------------------------------
  
    
  G4RotationMatrix rm3;
  rm3.rotateZ(180 *deg);
  
  G4ThreeVector positionMod3 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod3 = new G4Tubs("Mod3",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0/2., 
			  spanningAngleOfTheTube0);
  
  logicMod3 = new G4LogicalVolume(solidMod3, 
                                  Mod0Mater, "Mod3",0,0,0);
  
  
  physiMod3 = new G4PVPlacement(G4Transform3D(rm3, positionMod3), 
				 logicMod3,    // its logical volume			  
				 "Mod3",        // its name
				 logicMotherMod,  // its mother  volume
				 false,         // no boolean operations
				 0);            // no particular field
  



 for (G4int i=1;i<StepNumbers;i++)
  {

  physiMod[2*(StepNumbers-1)+i] = new G4PVPlacement(0,               
				PositionMod[i],  
				logicMod[i],     
				"Modstep",        
				logicMod3,     
				false,         
				0);  
			
  }

  //----------------------------------------------------------
  // Mother volume of the fourth modulator quarter
  //----------------------------------------------------------
  
    
  G4RotationMatrix rm4;
  rm4.rotateZ(270 *deg);
  
  G4ThreeVector positionMod4 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod4 = new G4Tubs("Mod4",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0/2., 
			  spanningAngleOfTheTube0);
  
  logicMod4 = new G4LogicalVolume(solidMod4, 
                                  Mod0Mater, "Mod4",0,0,0);
  
  
  physiMod4 = new G4PVPlacement(G4Transform3D(rm4, positionMod4), 
				 logicMod4,    			  
				 "Mod4",        
				 logicMotherMod,  
				 false,         
				 0);           
  

for (G4int i=1;i<StepNumbers;i++)
  {
  physiMod[3*(StepNumbers-1)+i] = new G4PVPlacement(0,               
				PositionMod[i],  
				logicMod[i],     
				"Modstep",        
				logicMod4,     
				false,         
				0);   
  }
  // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
  G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  red-> SetForceSolid(true);
  logicMotherMod -> SetVisAttributes(G4VisAttributes::Invisible);
 
  logicMod1 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod2 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod3 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod4 ->SetVisAttributes(G4VisAttributes::Invisible);
  
  for (G4int i=1;i<StepNumbers;i++)
  {
	 logicMod[i] -> SetVisAttributes(red); 
  }
  
 }

/////////////////////////////////////////////////////////////////////////////
// Messenger values
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator::SetModulatorAngle(G4double angle)
{
  G4double rotationAngle = angle;
  rm -> rotateZ(rotationAngle);
  physiMotherMod -> SetRotation(rm);  
  G4cout << "MODULATOR HAS BEEN ROTATED OF " << rotationAngle/deg 
	 << " deg" << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
}
/////////////////////////////////////////////////////////////////////////
// Change modulator material
void HadrontherapyModulator::SetModulatorMaterial(G4String Material)
{
    if (G4Material* NewMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Material, false) )
    {
	if (NewMaterial) 
	{
	    for(G4int i=1;i<StepNumbers;i++)
	    {
	    logicMod[i] -> SetMaterial(NewMaterial);  
	  //  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    
	  //  G4cout<<(logicMod[i]->GetMaterial()->GetName())<<G4endl;
	}
	G4cout << "The material of the Modulator wheel has been changed to " << Material << G4endl;
	}
    }
    else
    {
	G4cout << "WARNING: material \"" << Material << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl; 
	
	
    }
}

////////////////////////////////////////////////////////////////////////////////
// Change modulator position in the beam line
void HadrontherapyModulator::SetModulatorPosition(G4ThreeVector Pos)
{
  G4ThreeVector NewModulatorPos=Pos; 
  physiMotherMod -> SetTranslation( NewModulatorPos); 
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "The modulator wheel is translated to"<<  NewModulatorPos/mm <<"mm " <<G4endl;
	
}
/////////////////////////////////////////////////////////////////////////////////
//change modulator inner raduis
void HadrontherapyModulator::SetModulatorInnerRadius(G4double newvalue)
{
solidMod1 -> SetInnerRadius(newvalue);
solidMod2 -> SetInnerRadius(newvalue);
solidMod3 -> SetInnerRadius(newvalue);
solidMod4 -> SetInnerRadius(newvalue);
 for(G4int i=1;i<StepNumbers;i++)
 {
 solidMod[i] -> SetInnerRadius(newvalue);}
   G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
  G4cout << "InnerRadius of the Modulator Wheel has been changed to :"
	 << newvalue/mm<<" mm"<< G4endl;	 	
}
/////////////////////////////////////////////////////////////////////////////////
//change modulator outer raduis
void HadrontherapyModulator::SetModulatorOuterRadius(G4double newvalue)
{
solidMod1 -> SetOuterRadius(newvalue);
solidMod2 -> SetOuterRadius(newvalue);
solidMod3 -> SetOuterRadius(newvalue);
solidMod4 -> SetOuterRadius(newvalue);
 for(G4int i=1;i<StepNumbers;i++)
 {
 solidMod[i] -> SetOuterRadius(newvalue);}
   G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
  G4cout << "OuterRadius of the Modulator Wheel has been changed to :"
	 << newvalue/mm<<" mm"<<G4endl;	
}
/////////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator:: GetDataFromFile(G4String value)

{
G4String Name=value;
if(value=="default" )	
{
Name=FileName;
}
G4cout<<" Step properties of modulator will be get out from the external file "
 <<Name<<G4endl;
ModulatorPropertiesFromFile(Name);
}
