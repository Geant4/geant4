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
/// \file B107FastSim/src/G4ChannelingFastSimModel.cc
/// \brief Implementation of the G4ChannelingFastSimModel class
//
//
//
#include "G4ChannelingFastSimModel.hh"

#include "Randomize.hh"

#include "G4TransportationManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingFastSimModel::G4ChannelingFastSimModel(const G4String& modelName, G4Region* envelope)
: G4VFastSimulationModel(modelName, envelope)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingFastSimModel::G4ChannelingFastSimModel(const G4String& modelName)
: G4VFastSimulationModel(modelName)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingFastSimModel::~G4ChannelingFastSimModel()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4ChannelingFastSimModel::IsApplicable(const G4ParticleDefinition& particleType)
{
  return std::abs(particleType.GetPDGCharge())>DBL_EPSILON;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4ChannelingFastSimModel::ModelTrigger(const G4FastTrack& fastTrack)
{
  //default output
  G4bool modelTrigger = false;

  G4int particleDefinitionID =
          fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleDefinitionID();
  //kinetic energy
  G4double ekinetic = fastTrack.GetPrimaryTrack()->GetKineticEnergy();

  //energy cut, at the beginning, to not check everything else
  if(ekinetic > GetLowKineticEnergyLimit(particleDefinitionID))
  {
      //current logical volume
      G4LogicalVolume* crystallogic = fastTrack.GetEnvelopeLogicalVolume();
      fCrystalData->SetGeometryParameters(crystallogic);

      G4ThreeVector momentumDirection = fastTrack.GetPrimaryTrackLocalDirection();
      // the particle angle vs crystal plane or axis
      G4double angle = std::atan(momentumDirection.x()/momentumDirection.z());
      //recalculate angle into the lattice reference system
      angle = fCrystalData->
              AngleXFromBoxToLattice(angle,
                                     (fCrystalData->CoordinatesFromBoxToLattice(
                                          fastTrack.GetPrimaryTrackLocalPosition())).z());
      if (fCrystalData->GetModel()==2)
      {
          angle = std::sqrt(angle*angle+
                            std::pow(std::atan(momentumDirection.y()/
                                               momentumDirection.z()),2));
      }

      //particle mass
      G4double mass = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGMass();
      //particle total energy
      G4double etotal = fastTrack.GetPrimaryTrack()->GetTotalEnergy();

      //Particle position
      G4ThreeVector xyz0 = fastTrack.GetPrimaryTrackLocalPosition();
      //Step estimate
      G4double dz0 = fCrystalData->GetMaxSimulationStep(etotal,mass);
      xyz0 += 2*dz0*momentumDirection;//overestimated particle shift on the next step
                                      //in channeling

      //Applies the parameterisation not at the last step, only forward local direction
      //above low energy limit and below angular limit

      modelTrigger = (crystallogic->GetSolid()->
                      Inside(xyz0)==kInside) &&
                      momentumDirection.z()>0. &&
                      std::abs(angle) < GetLindhardAngleNumberHighLimit(particleDefinitionID) *
                      fCrystalData->GetLindhardAngle(etotal,mass);
  }

  return modelTrigger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChannelingFastSimModel::DoIt(const G4FastTrack& fastTrack,
                     G4FastStep& fastStep)
{  
  G4double etotal;//particle total energy
  G4double etotalPreStep;//etotal at the previous step
  G4double etotalToSetParticleProperties;//etotal value at which
                                         //SetParticleProperties is called
  G4double mass;  //particle mass
  G4double charge;//particle charge
  G4double tGlobal; //global time
  G4double tGlobalPreStep; //global time at the previous step
  G4ThreeVector xyz0;// the coordinates in the local reference system of the volume
  G4ThreeVector xyz0PreStep;// xyz at the previous step
  G4ThreeVector xyz;// the coordinates in the co-rotating reference system within
                    //a channel (elementary periodic cell)
  G4double x,y,z;   // the coordinates in the co-rotating reference system within
                    //a channel (elementary periodic cell)
  G4double tx0,ty0; // the angles in the local reference system of the volume
  G4double tx,ty;   // the angles in the co-rotating reference system within
                    //a channel (elementary periodic cell)
  G4double txPreStep,tyPreStep;// tx,ty at the previous step
  G4ThreeVector momentumDirection;
  G4ThreeVector scatteringAnglesAndEnergyLoss;//output of scattering functions
  G4double lindhardAngleNumberHighLimit0; //current high limit of the angle expressed in
                                          //[Lindhard angle] units

  //coordinates in Runge-Kutta calculations
  G4double x1=0.,x2=0.,x3=0.,x4=0.,y1=0.,y2=0.,y3=0.,y4=0.;
  //angles in Runge-Kutta calculations
  G4double tx1=0.,tx2=0.,tx3=0.,tx4=0.,ty1=0.,ty2=0.,ty3=0.,ty4=0.;
  //variables in Runge-Kutta calculations
  G4double kvx1=0.,kvx2=0.,kvx3=0.,kvx4=0.,kvy1=0.,kvy2=0.,kvy3=0.,kvy4=0.;
  //simulation step along z (internal step of the model) and its parts
  G4double dz,dzd3,dzd8;//dzd3 = dz/3; dzd8 = dz/8;
  //simulation step along the momentum direction
  G4double momentumDirectionStep;
  //effective simulation step (taking into account nuclear density along the trajectory)
  G4double effectiveStep;

  // flag, if Inside(xyz0) switches to kInside
  G4bool inside = false;

  G4LogicalVolume* crystallogic = fastTrack.GetEnvelopeLogicalVolume();
  fCrystalData->SetGeometryParameters(crystallogic);

  //set the max number of secondaries (photons) that can be added at this fastStep
  if (fRad)
  {
      fastStep.SetNumberOfSecondaryTracks(fMaxPhotonsProducedPerStep);
      //reseting the BaierKatkov integral to start it with the new trajectory
      fBaierKatkov->ResetRadIntegral();//to avoid any memory from the previous trajectory
  }

  etotal = fastTrack.GetPrimaryTrack()->GetTotalEnergy();
  mass = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGMass();
  charge = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGCharge();

  // we need to distunguish only charge particles, either leptons or hadrons
  G4bool hadron =
          fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetLeptonNumber()==0;

  lindhardAngleNumberHighLimit0 =
      GetLindhardAngleNumberHighLimit(fastTrack.GetPrimaryTrack()->
                                      GetParticleDefinition()->GetParticleDefinitionID());

  //set fCrystalData parameters depending on the particle parameters
  fCrystalData->SetParticleProperties(etotal, mass, charge, hadron);

  //global time
  tGlobal = fastTrack.GetPrimaryTrack()->GetGlobalTime();

  //coordinates in the co-rotating reference system within a channel
  xyz0= fastTrack.GetPrimaryTrackLocalPosition();
  xyz = fCrystalData->CoordinatesFromBoxToLattice(xyz0);
  x=xyz.x();
  y=xyz.y();
  z=xyz.z();

  momentumDirection=fastTrack.GetPrimaryTrackLocalDirection();
  //angle in the co-rotating reference system within a channel
  //(!!! ONLY FORWARD DIRECTION, momentumDirection.getZ()>0,
  //valid for high energies defined by the standard energy cuts)
  tx0 = std::atan(momentumDirection.x()/momentumDirection.z());
  ty0 = std::atan(momentumDirection.y()/momentumDirection.z());

  //angles in the co-rotating reference system within a channel
  tx = fCrystalData->AngleXFromBoxToLattice(tx0,z);
  ty = ty0;

  etotalToSetParticleProperties = etotal*0.999;
  G4bool inCrystal=true;//flag necessary to escape the cycle (at inCrystal=0;)
  //do calculations until the particle is inside the volume
  do
  {
      //remember the global time before the next step dz
      tGlobalPreStep=tGlobal;
      //remember the coordinates before the next step dz
      xyz0PreStep = xyz0;
      //remember the angles and the total energy before the step dz
      txPreStep = tx;
      tyPreStep = ty;
      etotalPreStep = etotal;

      dz = fCrystalData->GetSimulationStep(tx,ty);
      dzd3=dz/3;
      dzd8=dz/8;

      //trajectory calculation:
      //Runge-Cutt "3/8"
      //fCrystalData->GetCurv()*fCrystalData->GetCorrectionZ() is due to dependence of
      //the radius on x; GetCurv gets 1/R for the central ("central plane/axis")

      //first step
      kvx1=fCrystalData->Ex(x,y);
      x1=x+tx*dzd3;
      tx1=tx+(kvx1-fCrystalData->GetCurv()*fCrystalData->GetCorrectionZ())*dzd3;
      if (fCrystalData->GetModel()==2)
      {
         kvy1=fCrystalData->Ey(x,y);
         y1=y+ty*dzd3;
         ty1=ty+kvy1*dzd3;
      }

      //second step
      kvx2=fCrystalData->Ex(x1,y1);
      x2=x-tx*dzd3+tx1*dz;
      tx2=tx-(kvx1-fCrystalData->GetCurv()*fCrystalData->GetCorrectionZ())*dzd3+
              (kvx2-fCrystalData->GetCurv()*fCrystalData->GetCorrectionZ())*dz;
      if (fCrystalData->GetModel()==2)
      {
         kvy2=fCrystalData->Ey(x1,y1);
         y2=y-ty*dzd3+ty1*dz;
         ty2=ty-kvy1*dzd3+kvy2*dz;
      }

      //third step
      kvx3=fCrystalData->Ex(x2,y2);
      x3=x+(tx-tx1+tx2)*dz;
      tx3=tx+(kvx1-kvx2+kvx3-fCrystalData->GetCurv()*fCrystalData->GetCorrectionZ())*dz;
      if (fCrystalData->GetModel()==2)
      {
         kvy3=fCrystalData->Ey(x2,y2);
         y3=y+(ty-ty1+ty2)*dz;
         ty3=ty+(kvy1-kvy2+kvy3)*dz;
      }

      //fourth step
      kvx4=fCrystalData->Ex(x3,y3);
      x4=x+(tx+3.*tx1+3.*tx2+tx3)*dzd8;
      tx4=tx+(kvx1+3.*kvx2+3.*kvx3+kvx4)*dzd8-
              fCrystalData->GetCurv()*fCrystalData->GetCorrectionZ()*dz;
      if (fCrystalData->GetModel()==2)
      {
          kvy4=fCrystalData->Ey(x3,y3);
          y4=y+(ty+3.*ty1+3.*ty2+ty3)*dzd8;
          ty4=ty+(kvy1+3.*kvy2+3.*kvy3+kvy4)*dzd8;
      }
      else
      {
          y4 =y+ty*dz;
          ty4=ty;
      }

      x=x4;
      tx=tx4;
      y=y4;
      ty=ty4;

      z+=dz*fCrystalData->GetCorrectionZ();//motion along the z coordinate
                                          //("central plane/axis", no current plane/axis)

      xyz = fCrystalData->ChannelChange(x,y,z);
      x=xyz.x();
      y=xyz.y();
      z=xyz.z();

      //the coordinates in the local reference system of the volume
      //this vector will be used in the cycle escape condition and
      //in the radiation model (if activated)
      xyz0=fCrystalData->CoordinatesFromLatticeToBox(xyz);

      momentumDirectionStep=
              dz*std::sqrt(1+std::pow(std::tan(tx),2)+std::pow(std::tan(ty),2));
      tGlobal+=momentumDirectionStep/(fCrystalData->GetBeta())/CLHEP::c_light;

      //default scattering and energy loss 0
      scatteringAnglesAndEnergyLoss = G4ThreeVector(0.,0.,0.);

      //calculate separately for each element of the crystal
      for (G4int i = 0; i < fCrystalData->GetNelements(); i++)
      {
          //effective step taking into account nuclear density along the trajectory
          effectiveStep = momentumDirectionStep*fCrystalData->NuclearDensity(x,y,i);
          //Coulomb scattering on screened atomic potential (both multiple and single)
          scatteringAnglesAndEnergyLoss += fCrystalData->
                         CoulombAtomicScattering(effectiveStep,momentumDirectionStep,i);

          //Amorphous part of ionization energy losses
          etotal-=fCrystalData->IonizationLosses(momentumDirectionStep, i);
      }
      //electron scattering and coherent part of ionization energy losses
      scatteringAnglesAndEnergyLoss += fCrystalData->CoulombElectronScattering(
                                                   fCrystalData->MinIonizationEnergy(x,y),
                                                   fCrystalData->ElectronDensity(x,y),
                                                   momentumDirectionStep);
      tx += scatteringAnglesAndEnergyLoss.x();
      ty += scatteringAnglesAndEnergyLoss.y();
      etotal -= scatteringAnglesAndEnergyLoss.z();

      // recalculate the energy depended parameters
      //(only if the energy decreased enough, not at each step)
      if (etotalToSetParticleProperties>etotal)
      {
          fCrystalData->SetParticleProperties(etotal, mass, charge, hadron);
          etotalToSetParticleProperties = etotal*0.999;
      }

      //chain of conditions to escape the cycle
      // if Inside(xyz0)==kInside has been already true
      //(a particle has been inside the crystal)
      if (inside)
      {
          // if low energy
         if (etotal-mass<=GetLowKineticEnergyLimit(fastTrack.GetPrimaryTrack()->
                                                   GetParticleDefinition()->
                                                   GetParticleDefinitionID()))
             {inCrystal = false;}//escape the cycle
         //check if the angle w.r.t. the axes or planes is too high =>
         //return to standard Geant4:
         else if (fCrystalData->GetModel()==1) //1D model, field of planes
         {
             //if the angle w.r.t. the planes is too high
             if (std::abs(tx) >=
                     lindhardAngleNumberHighLimit0*fCrystalData->GetLindhardAngle())
                {inCrystal = false;}//escape the cycle
         }
         else if (fCrystalData->GetModel()==2) //2D model, field of axes
         {
             //if the angle w.r.t. the axes is too high
             if (std::sqrt(tx*tx+ty*ty) >= lindhardAngleNumberHighLimit0*
                                           fCrystalData->GetLindhardAngle())
                {inCrystal = false;}//escape the cycle
         }

           //radiation production & radiation energy losses
           //works only if the radiation model is activated
           if (fRad)
           {
               //back to the local reference system of the volume
               tx0 = fCrystalData->AngleXFromLatticeToBox(tx,z);
               ty0 = ty;
               //xyz0 was calculated above

               //running the radiation model and checking if a photon has been emitted
               if(fBaierKatkov->DoRadiation(etotal,mass,
                                           tx0,ty0,
                                           scatteringAnglesAndEnergyLoss.x(),
                                           scatteringAnglesAndEnergyLoss.y(),
                                           momentumDirectionStep,tGlobal,xyz0,
                                           crystallogic->
                                           GetSolid()->
                                           Inside(xyz0)!=kInside&&inCrystal))
               // also it was checked if the particle is escaping the volume
               // calculate the radiation integral immidiately in this case
               {
                   //a photon has been emitted!
                   //shift the particle back into the radiation point
                   etotal = fBaierKatkov->GetParticleNewTotalEnergy();
                   tx0 = fBaierKatkov->GetParticleNewAngleX();
                   ty0 = fBaierKatkov->GetParticleNewAngleY();
                   tGlobal = fBaierKatkov->GetNewGlobalTime();
                   xyz0 = fBaierKatkov->GetParticleNewCoordinateXYZ();

                   //add secondary photon
                   fBaierKatkov->GeneratePhoton(fastStep);

                   //particle energy was changed
                   fCrystalData->SetParticleProperties(etotal, mass, charge, hadron);

                   //coordinates in the co-rotating reference system within a channel
                   xyz = fCrystalData->CoordinatesFromBoxToLattice(xyz0);
                   x=xyz.x();
                   y=xyz.y();
                   z=xyz.z();

                   //angles in the co-rotating reference system within a channel
                   tx = fCrystalData->AngleXFromBoxToLattice(tx0,z);
                   ty = ty0;
               }
           }

           //precise check if the particle is escaping the volume
           if (crystallogic->GetSolid()->
                   Inside(xyz0)!=kInside)
           {
               //one step back to remain inside the volume
               //after the escape of the volume
               tGlobal = tGlobalPreStep;
               xyz0 = xyz0PreStep;
               tx = txPreStep;
               ty = tyPreStep;
               etotal = etotalPreStep;
               z-=dz*fCrystalData->GetCorrectionZ();
               // change the flag => this particle will not enter
               // the model before escape this volume

               inCrystal = false; //escape the cycle
           }
      }
      else
      {
          // if Inside(xyz0)==kInside we can enable checking of particle escape
          if (crystallogic->GetSolid()->
                           Inside(xyz0)==kInside)
                 {inside = true;}
          // a very rare case, if a particle remains
          // on the boundary and escapes the crystal
          else if (crystallogic->GetSolid()->
                           Inside(xyz0)==kOutside)
                 {inCrystal = false;}//escape the cycle
      }
  }
  while (inCrystal);

  //the angles in the local reference system of the volume
  tx0 = fCrystalData->AngleXFromLatticeToBox(tx,z);
  ty0 = ty;

  //set global time
  fastStep.ProposePrimaryTrackFinalTime(tGlobal);
  //set final position
  fastStep.ProposePrimaryTrackFinalPosition(xyz0);
  //set final kinetic energy
  fastStep.ProposePrimaryTrackFinalKineticEnergy(etotal-
                         fastTrack.GetPrimaryTrack()->
                                   GetParticleDefinition()->GetPDGMass());
  //set final momentum direction
  G4double momentumDirectionZ =
          1./std::sqrt(1.+std::pow(std::tan(tx0),2)+std::pow(std::tan(ty0),2));
  fastStep.ProposePrimaryTrackFinalMomentumDirection(
              G4ThreeVector(momentumDirectionZ*std::tan(tx0),
                            momentumDirectionZ*std::tan(ty0),
                            momentumDirectionZ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChannelingFastSimModel::Input(const G4Material *crystal, const G4String &lattice)
{
   //initializing the class with containing all
   //the crystal material and crystal lattice data and
   //Channeling scattering and ionization processes
   fCrystalData = new G4ChannelingFastSimCrystalData();
   //setting all the crystal material and lattice data
   fCrystalData->SetMaterialProperties(crystal,lattice);

   //setting default low energy cuts for kinetic energy
   SetLowKineticEnergyLimit(1*GeV,"proton");
   SetLowKineticEnergyLimit(1*GeV,"anti_proton");
   SetLowKineticEnergyLimit(200*MeV,"e-");
   SetLowKineticEnergyLimit(200*MeV,"e+");

   //set the model high limit of the angle expressed in [Lindhard angle] units
   SetLindhardAngleNumberHighLimit(100.,"proton");
   SetLindhardAngleNumberHighLimit(100.,"anti_proton");
   SetLindhardAngleNumberHighLimit(100.,"e-");
   SetLindhardAngleNumberHighLimit(100.,"e+");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChannelingFastSimModel::RadiationModelActivate()
{
    fRad = true;
    //activate the Baier-Katkov radiation model
    fBaierKatkov = new G4BaierKatkov();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
