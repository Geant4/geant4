// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VClusterModel.cc,v 1.2 2000-07-20 12:02:38 grichine Exp $
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
  for( G4int i = 0 ; i < fClusterEnergyVector.entries() ; i++ )
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




