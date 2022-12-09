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
//
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ---------------- GFlashShowerModel ----------------
//
// Authors: E.Barberio & Joanna Weng - 9.11.2004
// ------------------------------------------------------------

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4PionZero.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "geomdefs.hh"

#include "GFlashShowerModel.hh"
#include "GFlashHomoShowerParameterisation.hh"
#include "GFlashSamplingShowerParameterisation.hh"
#include "GFlashEnergySpot.hh"


GFlashShowerModel::GFlashShowerModel(G4String modelName,
                                     G4Envelope* envelope)
  : G4VFastSimulationModel(modelName, envelope),
    PBound(0), Parameterisation(0), HMaker(0)
{
  FlagParamType           = 0;
  FlagParticleContainment = 1;  
  StepInX0 = 0.1;
  Messenger       = new GFlashShowerModelMessenger(this);
}

GFlashShowerModel::GFlashShowerModel(G4String modelName)
  : G4VFastSimulationModel(modelName),
    PBound(0), Parameterisation(0), HMaker(0)
{
  FlagParamType           =1;
  FlagParticleContainment = 1;  
  StepInX0 = 0.1; 
  Messenger       = new GFlashShowerModelMessenger(this); 
}

GFlashShowerModel::~GFlashShowerModel()
{
  delete Messenger;
}

G4bool
GFlashShowerModel::IsApplicable(const G4ParticleDefinition& particleType)
{ 
  return 
  &particleType == G4Electron::ElectronDefinition() ||
  &particleType == G4Positron::PositronDefinition(); 
}

/**********************************************************************/
/* Checks whether conditions of fast parameterisation  are fullfilled */
/**********************************************************************/

G4bool GFlashShowerModel::ModelTrigger(const G4FastTrack & fastTrack )

{
  G4bool select = false;
  if(FlagParamType != 0)                  
  {
    G4double  ParticleEnergy = fastTrack.GetPrimaryTrack()->GetKineticEnergy(); 
    G4ParticleDefinition &ParticleType =
      *(fastTrack.GetPrimaryTrack()->GetDefinition()); 
    if(ParticleEnergy > PBound->GetMinEneToParametrise(ParticleType) &&
       ParticleEnergy < PBound->GetMaxEneToParametrise(ParticleType) )
    {
      // check conditions depending on particle flavour
      // performance to be optimized @@@@@@@
      Parameterisation->GenerateLongitudinalProfile(ParticleEnergy);
      select     = CheckParticleDefAndContainment(fastTrack);  
      if (select) EnergyStop= PBound->GetEneToKill(ParticleType);
    }
  }

  return select; 
}


G4bool
GFlashShowerModel::CheckParticleDefAndContainment(const G4FastTrack& fastTrack)
{  
  G4bool filter=false;
  G4ParticleDefinition * ParticleType =
    fastTrack.GetPrimaryTrack()->GetDefinition(); 
  
  if(  ParticleType == G4Electron::ElectronDefinition() || 
    ParticleType == G4Positron::PositronDefinition() )
  {
    filter=true;
    if(FlagParticleContainment == 1)  
    {
      filter=CheckContainment(fastTrack); 
    }
  }
  return filter;  
}

G4bool GFlashShowerModel::CheckContainment(const G4FastTrack& fastTrack)
{
  G4bool filter=false;
  // track informations
  G4ThreeVector DirectionShower=fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector InitialPositionShower=fastTrack.GetPrimaryTrackLocalPosition();

  G4ThreeVector OrthoShower, CrossShower; 
  // Returns orthogonal vector 
  OrthoShower = DirectionShower.orthogonal();
  // Shower in direction perpendicular to OrthoShower and DirectionShower
  CrossShower = DirectionShower.cross(OrthoShower);
  
  G4double  R     = Parameterisation->GetAveR90();
  G4double  Z     = Parameterisation->GetAveT90();
  G4int CosPhi[4] = {1,0,-1,0};
  G4int SinPhi[4] = {0,1,0,-1};
  
  G4ThreeVector Position;
  G4int NlateralInside=0;
  // pointer to solid we're in
  G4VSolid *SolidCalo = fastTrack.GetEnvelopeSolid();
  for(int i=0; i<4 ;i++)
  {
    // polar coordinates
    Position = InitialPositionShower       + 
    Z*DirectionShower           +
    R*CosPhi[i]*OrthoShower     +
    R*SinPhi[i]*CrossShower     ;
    
    if(SolidCalo->Inside(Position) != kOutside) 
      NlateralInside++;
  }
  
  // choose to parameterise or flag when all inetc...
  if(NlateralInside==4) filter=true;
  // std::cout << " points =   " <<NlateralInside << std::endl;
  return filter;
}


void
GFlashShowerModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep)
{
  // parametrise electrons
  if(fastTrack.GetPrimaryTrack()->GetDefinition()
     == G4Electron::ElectronDefinition() || 
     fastTrack.GetPrimaryTrack()->GetDefinition()
     == G4Positron::PositronDefinition() ) 
  ElectronDoIt(fastTrack,fastStep);
}

void
GFlashShowerModel::ElectronDoIt(const G4FastTrack& fastTrack,
                                      G4FastStep&  fastStep)
{
  // std::cout<<"--- ElectronDoit --- "<<std::endl;
  
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);
  fastStep.SetTotalEnergyDeposited(fastTrack.GetPrimaryTrack()->
                                   GetKineticEnergy());
  
  //-----------------------------
  // Get track parameters 
  //-----------------------------  
  //E,vect{p} and t,vec(x)
  G4double Energy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  
  // axis of the shower, in global reference frame:
  G4ThreeVector DirectionShower =
    fastTrack.GetPrimaryTrack()->GetMomentumDirection();
  G4ThreeVector OrthoShower, CrossShower;
  OrthoShower = DirectionShower.orthogonal();
  CrossShower = DirectionShower.cross(OrthoShower);
  
  //--------------------------------
  ///Generate longitudinal profile
  //--------------------------------
  Parameterisation->GenerateLongitudinalProfile(Energy);
    // performance iteration @@@@@@@
  
  ///Initialisation of long. loop variables
  G4VSolid *SolidCalo = fastTrack.GetEnvelopeSolid();
  G4ThreeVector pos   = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector dir   = fastTrack.GetPrimaryTrackLocalDirection();
  G4double Bound      = SolidCalo->DistanceToOut(pos,dir); 
  
  G4double Dz       = 0.00;     
  G4double ZEndStep = 0.00;
  
  G4double EnergyNow        = Energy;
  G4double EneIntegral      = 0.00;   
  G4double LastEneIntegral  = 0.00;   
  G4double DEne             = 0.00;
  
  G4double NspIntegral      = 0.00;   
  G4double LastNspIntegral  = 0.00;   
  G4double DNsp             = 0.00;
  
  // starting point of the shower:
  G4ThreeVector PositionShower  = fastTrack.GetPrimaryTrack()->GetPosition();
  G4ThreeVector NewPositionShower    = PositionShower;   
  G4double      StepLenght           = 0.00;
  
  //--------------------------
  /// Begin Longitudinal Loop
  //-------------------------
  
  do
  {  
    //determine step size=min(1Xo,next boundary)
    G4double stepLength = StepInX0*Parameterisation->GetX0();
    if(Bound < stepLength)
    { 
      Dz    = Bound;
      Bound = 0.00;
    }
    else
    { 
      Dz    = stepLength;
      Bound = Bound-Dz;
    }
    ZEndStep=ZEndStep+Dz;
    
    // Determine Energy Release in Step
    if(EnergyNow > EnergyStop)
    {
      LastEneIntegral  = EneIntegral;
      EneIntegral      = Parameterisation->IntegrateEneLongitudinal(ZEndStep);
      DEne             = std::min( EnergyNow,
                                   (EneIntegral-LastEneIntegral)*Energy);
      LastNspIntegral  = NspIntegral;
      NspIntegral      = Parameterisation->IntegrateNspLongitudinal(ZEndStep);
      DNsp             = std::max(1., std::floor( (NspIntegral-LastNspIntegral)
                                                 *Parameterisation->GetNspot() ));
    }
    // end of the shower
    else
    {    
      DEne = EnergyNow;
      DNsp = std::max(1., std::floor( (1.- NspIntegral)
                                     *Parameterisation->GetNspot() ));
    } 
    EnergyNow  = EnergyNow - DEne;
    
    // Apply sampling fluctuation - only in sampling calorimeters
    //
    GFlashSamplingShowerParameterisation* sp =
      dynamic_cast<GFlashSamplingShowerParameterisation*>(Parameterisation);
    if (sp)
    {
      G4double DEneSampling = sp->ApplySampling(DEne,Energy);
      DEne = DEneSampling;
    }

    //move particle in the middle of the step
    StepLenght        = StepLenght + Dz/2.00;  
    NewPositionShower = NewPositionShower + 
    StepLenght*DirectionShower;
    StepLenght        = Dz/2.00;
    
    //generate spots & hits:
    for (G4int i = 0; i < DNsp; ++i)
    { 
      GFlashEnergySpot Spot;      
      
      //Spot energy: the same for all spots
      Spot.SetEnergy( DEne / DNsp );
      G4double PhiSpot = Parameterisation->GeneratePhi(); // phi of spot
      G4double RSpot   = Parameterisation                 // radius of spot
                         ->GenerateRadius(i,Energy,ZEndStep-Dz/2.);

      // check reference-> may be need to introduce rot matrix @@@
      // Position: equally spaced in z
      
      G4ThreeVector SpotPosition = NewPositionShower  +
            Dz/DNsp*DirectionShower*(i+1/2.-DNsp/2.)  +
            RSpot*std::cos(PhiSpot)*OrthoShower       +  
            RSpot*std::sin(PhiSpot)*CrossShower;      
      Spot.SetPosition(SpotPosition);
      
      //Generate Hits of this spot      
      HMaker->make(&Spot, &fastTrack);
    }
  }
  while(EnergyNow > 0.0 && Bound> 0.0);     
  
  //---------------
  /// End Loop
  //--------------- 
}

/*

void
GFlashShowerModel::GammaDoIt(const G4FastTrack& fastTrack,
                                   G4FastStep&  fastStep)
{ 
  
  if( fastTrack.GetPrimaryTrack()->GetKineticEnergy() > EnergyStop )
    return;
  
  //deposita in uno spot unico l'energia 
  //con andamento exp decrescente. 
  
  // Kill the particle to be parametrised
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);
  fastStep.SetTotalEnergyDeposited(fastTrack.GetPrimaryTrack()
                                   ->GetKineticEnergy());
  // other settings????
  feSpotList.clear(); 
  
  //-----------------------------
  // Get track parameters 
  //-----------------------------

  // E,vect{p} and t,vec(x)
  G4double Energy = 
    fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  // axis of the shower, in global reference frame:
  G4ThreeVector DirectionShower =
    fastTrack.GetPrimaryTrack()->GetMomentumDirection();
  // starting point of the shower:
  G4ThreeVector PositionShower =
    fastTrack.GetPrimaryTrack()->GetPosition();
  
  //G4double DEneSampling = Parameterisation->ApplySampling(Energy,Energy);
  //if(DEneSampling <= 0.00) DEneSampling=Energy;  
  
  if(Energy > 0.0)
  {
    G4double dist = Parameterisation->GenerateExponential(Energy);      
    
    GFlashEnergySpot Spot;
    Spot.SetEnergy( Energy );
    G4ThreeVector SpotPosition = PositionShower + dist*DirectionShower;  
    Spot.SetPosition(SpotPosition);
    
    // Record the Spot:
    feSpotList.push_back(Spot);
    
    //Generate Hits of this spot      
    HMaker->make(Spot);
  }
}

*/
