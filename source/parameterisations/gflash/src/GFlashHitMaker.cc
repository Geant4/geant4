// Created by  E.Barberio & Joanna Weng 9.11.2004

#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHandle.hh"

#include "GFlashHitMaker.hh"

GFlashHitMaker::GFlashHitMaker()
{
	fFakeStep          = new G4Step();
	fFakePreStepPoint  = fFakeStep->GetPreStepPoint();
	fFakePostStepPoint = fFakeStep->GetPostStepPoint();
	fTouchableHandle   = new G4TouchableHistory(); // talk to ?@@@
	fpNavigator        = new G4Navigator();
	fNaviSetup         = false;
}

GFlashHitMaker::~GFlashHitMaker()
{
	delete fFakeStep;
	delete fpNavigator;
}

void GFlashHitMaker::make(const GFlashEnergySpot &Spot)
{
	// Locate the spot
	if (!fNaviSetup)
	{
		fpNavigator->
		SetWorldVolume(G4TransportationManager::GetTransportationManager()->
		GetNavigatorForTracking()->GetWorldVolume() );
		fpNavigator->
		LocateGlobalPointAndUpdateTouchable(Spot.GetPosition(), fTouchableHandle(), false);
		fNaviSetup = true;
	}
	else
	{
		fpNavigator->
		LocateGlobalPointAndUpdateTouchable(Spot.GetPosition(), fTouchableHandle());
	}
	//--------------------------------------
	// Fills attribute of the G4Step needed
	// by our sensitive detector:
	//-------------------------------------
	// set spot information:
	fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
	fFakePreStepPoint->SetPosition(Spot.GetPosition());
	fFakeStep->SetTotalEnergyDeposit(Spot.GetEnergy());	
	//--------------------------------------
	// Produce Hits
	// call sensitive part: taken/adapted from the stepping:
	// Send G4Step information to Hit/Dig if the volume is sensitive
	//--------------G4TouchableHistory----------------------------------------
	
	G4VPhysicalVolume* pCurrentVolume = fFakeStep->GetPreStepPoint()->GetPhysicalVolume();		
	G4VSensitiveDetector* pSensitive;
	if( pCurrentVolume != 0 )
	{
		pSensitive = pCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
		if( pSensitive != 0 )
		{
			pSensitive->Hit(fFakeStep);
		}
		else // below for conservative programming
		{  
			#ifdef GFLASH_DEBUG
			std::cout << "Out of SD = "<<fFakeStep->GetPreStepPoint()->GetPosition()<<endl; 
			#endif
		}
	}
	else
	{     
		#ifdef GFLASH_DEBUG
		std::cout << "GFlashHitMaker::Out of volume  "<<endl;
		#endif
	}
}


