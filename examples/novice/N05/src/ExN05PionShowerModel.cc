// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05PionShowerModel.cc,v 1.1 1999-01-07 16:06:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "ExN05PionShowerModel.hh"
#include "ExN05EnergySpot.hh"

#include "Randomize.hh"

#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"

#include "G4Colour.hh"

ExN05PionShowerModel::ExN05PionShowerModel(G4String modelName, G4LogicalVolume* envelope)
: G4VFastSimulationModel(modelName, envelope)
{
  fFakeStep          = new G4Step();
  fFakePreStepPoint  = fFakeStep->GetPreStepPoint();
  fFakePostStepPoint = fFakeStep->GetPostStepPoint();
  fpTouchable = new G4TouchableHistory();
  fpNavigator = new G4Navigator();
  fNaviSetup  = false;
}

ExN05PionShowerModel::ExN05PionShowerModel(G4String modelName)
: G4VFastSimulationModel(modelName)
{
  fFakeStep          = new G4Step();
  fFakePreStepPoint  = fFakeStep->GetPreStepPoint();
  fFakePostStepPoint = fFakeStep->GetPostStepPoint();
  fpTouchable = new G4TouchableHistory();
  fpNavigator = new G4Navigator();
  fNaviSetup  = false;
}

ExN05PionShowerModel::~ExN05PionShowerModel()
{
  delete fFakeStep;
  delete fpTouchable;
  delete fpNavigator;
}

G4bool ExN05PionShowerModel::IsApplicable(const G4ParticleDefinition& particleType)
{
  return 
    &particleType == G4PionMinus::PionMinusDefinition() ||
    &particleType == G4PionPlus::PionPlusDefinition();
}

G4bool ExN05PionShowerModel::ModelTrigger(const G4FastTrack& fastTrack)
{
  // Applies the parameterisation always: ie as soon as the pion enters the envelope
  return true;
}

void ExN05PionShowerModel::DoIt(const G4FastTrack& fastTrack, 
		     G4FastStep& fastStep)
{
  // Kill the parawmeterised particle:
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);
  fastStep.SetTotalEnergyDeposited(fastTrack.GetPrimaryTrack()->GetKineticEnergy());

  // split into "energy spots" energy according to the shower shape:
  Explode(fastTrack);
  
  // and put those energy spots into the crystals:
  BuildDetectorResponse();
  
}

void ExN05PionShowerModel::Explode(const G4FastTrack& fastTrack)
{
  //-----------------------------------------------------
  // Non-physical shower generated: exp along z and
  // transverse.
  //-----------------------------------------------------

  // center of the shower, we put at the middle of the ghost:
  G4ThreeVector showerCenter;
  G4double distOut;
  distOut = fastTrack.GetEnvelopeSolid()->
    DistanceToOut(fastTrack.GetPrimaryTrackLocalPosition(),
		  fastTrack.GetPrimaryTrackLocalDirection());
  showerCenter = fastTrack.GetPrimaryTrackLocalPosition() + 
    (distOut/2.)*fastTrack.GetPrimaryTrackLocalDirection();

  showerCenter = fastTrack.GetInverseAffineTransformation()->
    TransformPoint(showerCenter);

  // axis of the shower, in global reference frame:
  G4ThreeVector xShower, yShower, zShower;
  zShower = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
  xShower = zShower.orthogonal();
  yShower = zShower.cross(xShower);
  
  // shoot the energy spots:
  G4double Energy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  G4int nSpot = 50;
  G4double deposit = Energy/double(nSpot);;
  ExN05EnergySpot eSpot;
  eSpot.SetEnergy(deposit);
  G4ThreeVector ePoint;

  G4double z, r, phi;
  for (int i = 0; i < nSpot; i++)
    {
      z   = RandGauss::shoot(0,20*cm);
      r   = RandGauss::shoot(0,10*cm);
      phi = RandFlat::shoot()*twopi;
      ePoint = showerCenter +
	z*zShower +
	r*cos(phi)*xShower + r*sin(phi)*yShower;
      eSpot.SetPosition(ePoint);
      feSpotList.insert(eSpot);
    }
}


void ExN05PionShowerModel::BuildDetectorResponse()
{
  // Does the assignation of the energy spots to the sensitive volumes:
  for (int i = 0; i < feSpotList.entries(); i++)
    {
      // Draw the energy spot:
      G4Colour red(1.,0.,0.);
      feSpotList[i].Draw(&red);
      //      feSpotList[i].Print();
      
      // "converts" the energy spot into the fake
      // G4Step to pass to sensitive detector:
      AssignSpotAndCallHit(feSpotList[i]);
    }
}


void ExN05PionShowerModel::AssignSpotAndCallHit(const ExN05EnergySpot &eSpot)
{
  //
  // "converts" the energy spot into the fake
  // G4Step to pass to sensitive detector:
  //
  FillFakeStep(eSpot);

  //
  // call sensitive part: taken/adapted from the stepping:
  // Send G4Step information to Hit/Dig if the volume is sensitive
  //
  G4VPhysicalVolume* pCurrentVolume = 
    fFakeStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VSensitiveDetector* pSensitive;
  
  if( pCurrentVolume != 0 )
    {
      pSensitive = pCurrentVolume->GetLogicalVolume()->
	GetSensitiveDetector();
      if( pSensitive != 0 )
	{
	  pSensitive->Hit(fFakeStep);
	}
    }
}


void ExN05PionShowerModel::FillFakeStep(const ExN05EnergySpot &eSpot)
{
  //-----------------------------------------------------------
  // find in which volume the spot is.
  //-----------------------------------------------------------
  if (!fNaviSetup)
    {
      fpNavigator->
	SetWorldVolume(G4TransportationManager::GetTransportationManager()->
		       GetNavigatorForTracking()->GetWorldVolume());
      fpNavigator->
	LocateGlobalPointAndUpdateTouchable(eSpot.GetPosition(),
					    fpTouchable,
					    false);
      fNaviSetup = true;
    }
  else
    {
      fpNavigator->
	LocateGlobalPointAndUpdateTouchable(eSpot.GetPosition(),
					    fpTouchable);
     }
  //--------------------------------------
  // Fills attribute of the G4Step needed
  // by our sensitive detector:
  //-------------------------------------
  // set touchable volume at PreStepPoint:
  fFakePreStepPoint->SetTouchable(fpTouchable);
  // set total energy deposit:
  fFakeStep->SetTotalEnergyDeposit(eSpot.GetEnergy());
}






