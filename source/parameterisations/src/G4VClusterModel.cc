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
// $Id: G4VClusterModel.cc,v 1.3 2001-09-18 09:02:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Timer.hh"

#include "G4VClusterModel.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"


////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4VClusterModel::G4VClusterModel(const G4String& modelName,G4Envelope *anEnvelope) :
  G4VFastSimulationModel(modelName,anEnvelope)
{
  fFakeStep = new G4Step() ;
  fFakePreStepPoint = fFakeStep->GetPreStepPoint() ;
  fFakePostStepPoint = fFakeStep->GetPostStepPoint() ;

  fTouchable = new G4TouchableHistory ;
  fNavigator = new G4Navigator ;
  fNavigatorSetup = false ;

}

///////////////////////////////////////////////////////////////////////////

G4VClusterModel::~G4VClusterModel()
{
  delete fFakeStep ;
  delete fTouchable ;
  delete fNavigator ; 
}


//////////////////////////////////////////////////////////////////////////////
//
// 

void G4VClusterModel::BuildDetectorResponse()
{
  for( size_t i = 0 ; i < fClusterEnergyVector.entries() ; i++ )
  {
    AssignClusterHit(fClusterPositionVector[i],fClusterEnergyVector[i]) ;
  }
  fClusterPositionVector.clear() ;
  fClusterEnergyVector.clear() ;
}

//////////////////////////////////////////////////////////////////////////////
//
// 

void G4VClusterModel::AssignClusterHit(const G4ThreeVector& position, 
                                               G4double energy           )
{
  // "converts" the energy spot into the fake
  // G4Step to pass to sensitive detector:
  
  FillFakeStep(position,energy) ;

  // call sensitive part: taken/adapted from the stepping:
  // Send G4Step information to Hit/Dig if the volume is sensitive
  
  G4VPhysicalVolume* pCurrentVolume = fFakeStep->GetPreStepPoint()->
                                      GetPhysicalVolume() ;
  G4VSensitiveDetector* pSensitive ;
  
  if( pCurrentVolume != 0 )
  {
    pSensitive = pCurrentVolume->GetLogicalVolume()->GetSensitiveDetector() ;

    if( pSensitive != 0 )  pSensitive->Hit(fFakeStep) ;
  }
}

//////////////////////////////////////////////////////////////////////////////
//
// 

void G4VClusterModel::FillFakeStep(const G4ThreeVector& position, 
                                               G4double energy           )
{
  // find in which volume the spot is.

  if (!fNavigatorSetup)
  {
    fNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->
		       GetNavigatorForTracking()->GetWorldVolume()) ;

    fNavigator->LocateGlobalPointAndUpdateTouchable(position,fTouchable,false);
    fNavigatorSetup = true;
  }
  else
  {
    fNavigator->LocateGlobalPointAndUpdateTouchable(position,fTouchable);
  }
  // Fills attribute of the G4Step needed
  // by our sensitive detector:
  //
  // set touchable volume at PreStepPoint:

  fFakePreStepPoint->SetTouchable(fTouchable) ;

  // set total energy deposit:

  fFakeStep->SetTotalEnergyDeposit(energy) ;

  return ;
}


//
//
///////////////////////////////////////////////////////////////////////




