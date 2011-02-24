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
// Created by Joanna Weng 26.11.2004
#include <iostream>
// G4 Classes
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Colour.hh"

// User Classes
#include "ExGflashDetectorConstruction.hh"
#include "ExGflashSensitiveDetector.hh"
#include "ExGflashMaterialManager.hh"
//fast simulation
#include "GFlashHomoShowerParameterisation.hh"
#include "G4FastSimulationManager.hh"
#include "GFlashShowerModel.hh"
#include "GFlashHitMaker.hh"
#include "GFlashParticleBounds.hh"

using namespace std;

ExGflashDetectorConstruction::ExGflashDetectorConstruction()
:m_experimentalHall_log(0), 
m_calo_log(0), 
m_experimentalHall_phys(0), 
m_calo_phys(0)
{
	G4cout<<"ExGflashDetectorConstruction::Detector constructor"<<G4endl;    
	
	// Simplified `CMS-like` PbWO4 crystal calorimeter  
	m_calo_xside=31*cm;
	m_calo_yside=31*cm;
	m_calo_zside=24*cm; 
	
	// GlashStuff
	m_theParticleBounds  = new GFlashParticleBounds();		// Energy Cuts to kill particles
	m_theHMaker          = new GFlashHitMaker();     		// Makes the EnergieSpots	
}


ExGflashDetectorConstruction::~ExGflashDetectorConstruction()
{ 
//@@@ ExGflashDetectorConstruction::Soll ich alles dlete
	delete m_theParameterisation;
	delete m_theParticleBounds;
	delete m_theHMaker;
	delete m_theFastShowerModel;
}


G4VPhysicalVolume* ExGflashDetectorConstruction::Construct()
{
	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
	G4String mat= "PbWO4";     
	G4cout<<"Defining the materials"<<G4endl;
	ExGflashMaterialManager *matManager=ExGflashMaterialManager::GetMaterialManager();
	
	/*******************************
	* The Experimental Hall       *
	*******************************/
	m_experimentalHall_x=1000.*cm;
	m_experimentalHall_y=1000.*cm;
	m_experimentalHall_z=1000.*cm;
	
	m_experimentalHall_box = new G4Box("expHall_box",              // World Volume
	m_experimentalHall_x,           // x size
	m_experimentalHall_y,           // y size 
	m_experimentalHall_z);          // z size
	
	m_experimentalHall_log = new G4LogicalVolume(m_experimentalHall_box, //its solid
	matManager->getMaterial("Air"),	 	//its material
	"expHall_log", 					 	// its name 
	0,     								//opt: fieldManager
	0,     								//opt: SensitiveDetector 
	0);    								//opt: UserLimits                       
	m_experimentalHall_phys = new G4PVPlacement(0,    	            //no rotation       
	G4ThreeVector(),  	    			//at (0,0,0)
	"expHall",	            			//its name 
	m_experimentalHall_log,  			//its logical volume
	0,	                   				//its mother  volume
	false,                				//no boolean operation
	0);	            	 				//copy number
	

	//------------------------------ 
	// Calorimeter segments
	//------------------------------
	// Simplified `CMS-like` PbWO4 crystal calorimeter  

	m_NbOfCrystals = 10;  // this are the crystals PER ROW in this example 
	                      // cube of 10 x 10 crystals 
	                      // don't change it @the moment, since 
	                      // the readout in event action assumes this 
	                      // dimensions and is not automatically adapted
	                      // in this version of the example :-( 
	m_CrystalWidht = 3*cm;
	m_CrystalLenght= 24*cm;
	m_calo_xside=(m_CrystalWidht*m_NbOfCrystals)+1*cm;
	m_calo_yside=(m_CrystalWidht*m_NbOfCrystals)+1*cm;;
	m_calo_zside=m_CrystalLenght;
	
	G4Box *calo_box= new G4Box("CMS calorimeter",   	//its name
	m_calo_xside/2., 				 					//size						   
	m_calo_yside/2.,
	m_calo_zside/2.);
	m_calo_log = new G4LogicalVolume(calo_box,			//its solid
	matManager->getMaterial("Air"),			//its material
	"calo log",	                			//its name
	0,   									//opt: fieldManager                            
	0,   									//opt: SensitiveDetector 
	0);  									//opt: UserLimit
	
	G4double Xpos = 0.0;
	G4double Ypos = 0.0;
	G4double Zpos = 100.0*cm;
	
	m_calo_phys = new G4PVPlacement(0,	//no rotation
	G4ThreeVector(Xpos,Ypos,Zpos),		//at (0,0,0)
	m_calo_log,							//its logical volume				     
	"calorimeter",						//its name
	m_experimentalHall_log,  			//its mother  volume
	false,                  			//no boolean operation
	1);
			       					//Visibility
	for (int i=0; i<m_NbOfCrystals;i++)
	{
		
		for (int j=0; j<m_NbOfCrystals;j++)
		{	
			int n =  i*10+j;
			m_crystal[n]= new G4Box("Crystal",						//its name
			m_CrystalWidht/2,m_CrystalWidht/2,m_CrystalLenght/2);  	//size			
			m_crystal_log[n] = new G4LogicalVolume(m_crystal[n],	//its solid
			matManager->getMaterial(mat),							//its material
			"Crystal_log");											//its name
	      
			m_crystal_phys[n] = new G4PVPlacement(0,				//no rotation
			G4ThreeVector((i*m_CrystalWidht)-135,(j*m_CrystalWidht)-135,0 ),			//at (0,0,0)
			m_crystal_log[n],										//its logical volume				     
			"crystal",												//its name
			m_calo_log,               								//its mother volume
			false,                  								//no boolean operation
			1);														//Visibility
		}
	}	
	G4cout << "There are " << m_NbOfCrystals << " crystals per row in the calorimeter, so in total "<<
	m_NbOfCrystals*m_NbOfCrystals << " crystals" << G4endl;	
	G4cout << "The have widthof  " << m_CrystalWidht /cm << "  cm and a lenght of  " <<  m_CrystalLenght /cm
	<<" cm. The Material is "<< matManager->getMaterial(mat) << G4endl;
	
	// Sensitive Detector part
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	ExGflashSensitiveDetector* CaloSD=
	new ExGflashSensitiveDetector("Calorimeter",this);
	SDman->AddNewDetector(CaloSD);
 	
	m_experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);					   
	G4VisAttributes* CaloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	G4VisAttributes* CrystalVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	m_calo_log->SetVisAttributes(CaloVisAtt);
	for (int i=0; i<100;i++)
	{
		m_crystal_log[i]->SetVisAttributes(CrystalVisAtt);
		m_crystal_log[i]->SetSensitiveDetector(CaloSD);
	}
	// define the parameterisation region
	aRegion = new G4Region("crystals");
	m_calo_log->SetRegion(aRegion);
	aRegion->AddRootLogicalVolume(m_calo_log);
  
	
	/**********************************************
	* Initializing shower modell
	***********************************************/	
	G4cout<<"Shower parameterization"<<G4endl;
	m_theFastShowerModel =  new GFlashShowerModel("fastShowerModel",aRegion);
	m_theParameterisation = new GFlashHomoShowerParameterisation(matManager->getMaterial(mat));
	m_theFastShowerModel->SetParameterisation(*m_theParameterisation);
	m_theFastShowerModel->SetParticleBounds(*m_theParticleBounds) ;
	m_theFastShowerModel->SetHitMaker(*m_theHMaker);	 
	G4cout<<"end shower parameterization"<<G4endl;
	/**********************************************/
	
	return m_experimentalHall_phys;
}














