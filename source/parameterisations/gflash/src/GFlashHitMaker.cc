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
// Created by  E.Barberio & Joanna Weng 9.11.2004

#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHandle.hh"
#include "G4VGFlashSensitiveDetector.hh"

#include "GFlashHitMaker.hh"
#include "G4GFlashSpot.hh"

GFlashHitMaker::GFlashHitMaker()
{
	fTouchableHandle   = new G4TouchableHistory(); // talk to ?@@@
	fpNavigator        = new G4Navigator();
	fNaviSetup         = false;
}

GFlashHitMaker::~GFlashHitMaker()
{
	delete fpNavigator;
}

void GFlashHitMaker::make(GFlashEnergySpot * aSpot, const G4FastTrack * aT)
{
	// Locate the spot
	if (!fNaviSetup)
	{
		fpNavigator->
		SetWorldVolume(G4TransportationManager::GetTransportationManager()->
		GetNavigatorForTracking()->GetWorldVolume() );
		fpNavigator->
		LocateGlobalPointAndUpdateTouchable(aSpot->GetPosition(), fTouchableHandle(), false);
		 fNaviSetup = true;
	}
	else
	{
		fpNavigator->
		LocateGlobalPointAndUpdateTouchable(aSpot->GetPosition(), fTouchableHandle());
	}
	
	//--------------------------------------
	// Fills attribute of the G4Step needed
	// by our sensitive detector:
	//-------------------------------------
	// set spot information:
	G4GFlashSpot theSpot(aSpot, aT, fTouchableHandle);
	///Navigator
	//--------------------------------------
	// Produce Hits
	// call sensitive part: taken/adapted from the stepping:
	// Send G4Step information to Hit/Dig if the volume is sensitive
	//--------------G4TouchableHistory----------------------------------------
	
	G4VPhysicalVolume* pCurrentVolume = fTouchableHandle()->GetVolume();		
	G4VSensitiveDetector* pSensitive;
	if( pCurrentVolume != 0 )
	{
		pSensitive = pCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
		G4VGFlashSensitiveDetector * gflashSensitive = 
		               dynamic_cast<G4VGFlashSensitiveDetector * > (pSensitive);
		if( gflashSensitive )
		{
			gflashSensitive->Hit(&theSpot);
		}
		else if ( (!gflashSensitive ) && 
		         ( pSensitive ) && 
						 ( pCurrentVolume->GetLogicalVolume()->GetFastSimulationManager() )  
					 ) // Using gflash without implementing the 
							     // gflashSensitive detector interface -> not allowed!
		
		{  	
			std::cout<<"When using GFlash in geant4, please implement the "<<std::endl;
			std::cout<<"G4VGFlashSensitiveDetector interface in addition to the"<<std::endl;
			std::cout<<"G4VSensitiveDetector interface in the relevant sensitive detector."<<std::endl;
			G4Exception("GFlashHitMaker: G4VGFlashSensitiveDetector interface not implemented");
		}
	}
	else
	{     
		#ifdef GFLASH_DEBUG
		std::cout << "GFlashHitMaker::Out of volume  "<<endl;
		#endif
	}
}
