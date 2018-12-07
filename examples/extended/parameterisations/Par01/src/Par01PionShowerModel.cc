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
/// \file Par01/src/Par01PionShowerModel.cc
/// \brief Implementation of the Par01PionShowerModel class
//
//
//
#include "Par01PionShowerModel.hh"
#include "Par01EnergySpot.hh"

#include "Randomize.hh"

#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"

#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01PionShowerModel::Par01PionShowerModel(G4String modelName, G4Region* envelope)
: G4VFastSimulationModel(modelName, envelope)
{
  fFakeStep          = new G4Step();
  fFakePreStepPoint  = fFakeStep->GetPreStepPoint();
  fFakePostStepPoint = fFakeStep->GetPostStepPoint();
  fTouchableHandle   = new G4TouchableHistory();
  fpNavigator        = new G4Navigator();
  fNaviSetup         = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01PionShowerModel::Par01PionShowerModel(G4String modelName)
: G4VFastSimulationModel(modelName)
{
  fFakeStep          = new G4Step();
  fFakePreStepPoint  = fFakeStep->GetPreStepPoint();
  fFakePostStepPoint = fFakeStep->GetPostStepPoint();
  fTouchableHandle   = new G4TouchableHistory();
  fpNavigator        = new G4Navigator();
  fNaviSetup         = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01PionShowerModel::~Par01PionShowerModel()
{
  delete fFakeStep;
  delete fpNavigator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par01PionShowerModel::IsApplicable(const G4ParticleDefinition& particleType)
{
  return 
    &particleType == G4PionMinus::PionMinusDefinition() ||
    &particleType == G4PionPlus::PionPlusDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par01PionShowerModel::ModelTrigger(const G4FastTrack&)
{
  // Applies the parameterisation always:
  // ie as soon as the pion enters the envelope
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01PionShowerModel::DoIt(const G4FastTrack& fastTrack, 
                     G4FastStep& fastStep)
{
  //  G4cout << "Par01PionShowerModel::DoIt" << G4endl;

  // Kill the parameterised particle:
  fastStep.KillPrimaryTrack();
  fastStep.ProposePrimaryTrackPathLength(0.0);
  fastStep.ProposeTotalEnergyDeposited(fastTrack.GetPrimaryTrack()->GetKineticEnergy());

  // split into "energy spots" energy according to the shower shape:
  Explode(fastTrack);
  
  // and put those energy spots into the crystals:
  BuildDetectorResponse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01PionShowerModel::Explode(const G4FastTrack& fastTrack)
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
  G4double deposit = Energy/double(nSpot);
  Par01EnergySpot eSpot;
  eSpot.SetEnergy(deposit);
  G4ThreeVector ePoint;

  // clear the spot list before use
  feSpotList.clear();

  G4double z, r, phi;
  for (int i = 0; i < nSpot; i++)
    {
      z   = G4RandGauss::shoot(0,20*cm);
      r   = G4RandGauss::shoot(0,10*cm);
      phi = G4UniformRand()*twopi;
      ePoint = showerCenter +
        z*zShower +
        r*std::cos(phi)*xShower + r*std::sin(phi)*yShower;
      eSpot.SetPosition(ePoint);
      feSpotList.push_back(eSpot);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01PionShowerModel::BuildDetectorResponse()
{
  // Does the assignation of the energy spots to the sensitive volumes:
  for (size_t i = 0; i < feSpotList.size(); i++)
    {
      // Draw the energy spot:
      //      G4Colour red(1.,0.,0.);
      //      feSpotList[i].Draw(&red);
      //      feSpotList[i].Print();
      
      // "converts" the energy spot into the fake
      // G4Step to pass to sensitive detector:
      AssignSpotAndCallHit(feSpotList[i]);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01PionShowerModel::AssignSpotAndCallHit(const Par01EnergySpot &eSpot)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01PionShowerModel::FillFakeStep(const Par01EnergySpot &eSpot)
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
        LocateGlobalPointAndUpdateTouchableHandle(eSpot.GetPosition(),
                                                  G4ThreeVector(0.,0.,0.),
                                                  fTouchableHandle,
                                                  false);
      fNaviSetup = true;
    }
  else
    {
      fpNavigator->
        LocateGlobalPointAndUpdateTouchableHandle(eSpot.GetPosition(),
                                                  G4ThreeVector(0.,0.,0.),
                                                  fTouchableHandle);
    }
  //--------------------------------------
  // Fills attribute of the G4Step needed
  // by our sensitive detector:
  //-------------------------------------
  // set touchable volume at PreStepPoint:
  fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
  // set total energy deposit:
  fFakeStep->SetTotalEnergyDeposit(eSpot.GetEnergy());
}
