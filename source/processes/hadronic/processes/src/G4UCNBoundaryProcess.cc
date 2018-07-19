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
// $Id: G4UCNBoundaryPrcess.cc 69576 2013-05-08 13:48:13Z gcosmo $
//
///////////////////////////////////////////////////////////////////////
// UCN BoundaryProcess Class Implementation
///////////////////////////////////////////////////////////////////////
//
// File:        G4UCNBoundaryProcess.cc
// Description: Discrete Process -- Boundary Process of UCN
// Version:     1.0
// Created:     2014-05-12
// Author:      Peter Gumplinger
//              adopted from Geant4UCN by Peter Fierlinger (9.9.04) and
//              Marcin Kuzniak (21.4.06)
//              1/v energy dependent absorption cross section
//              inside materials
// Updated:     2007 Extensions for the microroughness model by Stefan Heule
//
// mail:        gum@triumf.ca
//
///////////////////////////////////////////////////////////////////////

#include "G4UCNProcessSubType.hh"


#include "G4UCNBoundaryProcess.hh"
#include "G4UCNBoundaryProcessMessenger.hh"

#include "G4GeometryTolerance.hh"

#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"

#include "G4UCNMaterialPropertiesTable.hh"

#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4VSensitiveDetector.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4UCNBoundaryProcess::G4UCNBoundaryProcess(const G4String& processName,
                                           G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{

  if (verboseLevel > 0) G4cout << GetProcessName() << " is created " << G4endl;

  SetProcessSubType(fUCNBoundary);

  theStatus = Undefined;

  fMessenger = new G4UCNBoundaryProcessMessenger(this);

  neV = 1.0e-9*eV;

  kCarTolerance = G4GeometryTolerance::GetInstance()
                  ->GetSurfaceTolerance();

  Material1 = NULL;
  Material2 = NULL;

  aMaterialPropertiesTable1 = NULL;
  aMaterialPropertiesTable2 = NULL;

  UseMicroRoughnessReflection = false;
  DoMicroRoughnessReflection  = false;

  nNoMPT = nNoMRT = nNoMRCondition = 0;
  nAbsorption = nEzero = nFlip = 0;
  aSpecularReflection = bSpecularReflection = 0;
  bLambertianReflection = 0;
  aMRDiffuseReflection = bMRDiffuseReflection = 0;
  nSnellTransmit = mSnellTransmit = 0;
  aMRDiffuseTransmit = 0;
  ftheta_o = fphi_o = 0;
}

G4UCNBoundaryProcess::~G4UCNBoundaryProcess()
{
   delete fMessenger;
} 

G4VParticleChange*
G4UCNBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeVelocity(aTrack.GetVelocity());

  // Get hyperStep from  G4ParallelWorldProcess
  //  NOTE: PostSetpDoIt of this process should be
  //        invoked after G4ParallelWorldProcess!

  const G4Step* pStep = &aStep;

  const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();

  if (hStep) pStep = hStep;

  G4bool isOnBoundary =
          (pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);

  if (isOnBoundary) {
     Material1 = pStep->GetPreStepPoint()->GetMaterial();
     Material2 = pStep->GetPostStepPoint()->GetMaterial();
  } else {
     theStatus = NotAtBoundary;
     if ( verboseLevel > 1 ) BoundaryProcessVerbose();
     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if (aTrack.GetStepLength()<=kCarTolerance/2) {
     theStatus = StepTooSmall;
     if ( verboseLevel > 0 ) BoundaryProcessVerbose();
     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  aMaterialPropertiesTable1 = (G4UCNMaterialPropertiesTable*)Material1->
                                                 GetMaterialPropertiesTable();
  aMaterialPropertiesTable2 = (G4UCNMaterialPropertiesTable*)Material2->
                                                 GetMaterialPropertiesTable();

  G4String volnam1 = pStep->GetPreStepPoint() ->GetPhysicalVolume()->GetName();
  G4String volnam2 = pStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();

  if (verboseLevel > 0) {
     G4cout << " UCN at Boundary! " << G4endl;
     G4cout << " vol1: " << volnam1 << ", vol2: " << volnam2 << G4endl;
     G4cout << " Ekin:     " << aTrack.GetKineticEnergy()/neV <<"neV"<< G4endl;
     G4cout << " MomDir:   " << aTrack.GetMomentumDirection() << G4endl;
  }

  if (Material1 == Material2) {
     theStatus = SameMaterial;
     if ( verboseLevel > 0 ) BoundaryProcessVerbose();
     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4ThreeVector theGlobalPoint = pStep->GetPostStepPoint()->GetPosition();

  G4bool valid;
  //  Use the new method for Exit Normal in global coordinates,
  //    which provides the normal more reliably.

  // ID of Navigator which limits step

  G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
          G4TransportationManager::GetTransportationManager()->
                                   GetActiveNavigatorsIterator();

  G4ThreeVector theGlobalNormal =
             (iNav[hNavId])->GetGlobalExitNormal(theGlobalPoint,&valid);

  if (valid) {
     theGlobalNormal = -theGlobalNormal;
  }
  else
  {
     G4ExceptionDescription ed;
     ed << " G4UCNBoundaryProcess/PostStepDoIt(): "
        << " The Navigator reports that it returned an invalid normal"
        << G4endl;
     G4Exception("G4UCNBoundaryProcess::PostStepDoIt", "UCNBoun01",
                 EventMustBeAborted,ed,
                 "Invalid Surface Normal - Geometry must return valid surface normal");
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

  G4ThreeVector OldMomentum = aParticle->GetMomentumDirection();

  if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4OPTICAL_DEBUG
     G4ExceptionDescription ed;
     ed << " G4UCNBoundaryProcess/PostStepDoIt(): "
        << " theGlobalNormal points in a wrong direction. "
        << G4endl;
     ed << "    The momentum of the photon arriving at interface (oldMomentum)"
        << " must exit the volume cross in the step. " << G4endl;
     ed << "  So it MUST have dot < 0 with the normal that Exits the new volume (globalNormal)." << G4endl;
     ed << "  >> The dot product of oldMomentum and global Normal is " << OldMomentum*theGlobalNormal << G4endl;
     ed << "     Old Momentum  (during step)     = " << OldMomentum << G4endl;
     ed << "     Global Normal (Exiting New Vol) = " << theGlobalNormal << G4endl;
     ed << G4endl;
     G4Exception("G4UCNBoundaryProcess::PostStepDoIt", "UCNBoun02",
                 EventMustBeAborted,  // Or JustWarning to see if it happens repeatedbly on one ray
                 ed,
                "Invalid Surface Normal - Geometry must return valid surface normal pointing in the right direction");
#else
     theGlobalNormal = -theGlobalNormal;
#endif
  }

  G4ThreeVector theNeutronMomentum = aTrack.GetMomentum();

  G4double theMomentumNormal = theNeutronMomentum*theGlobalNormal;
  G4double theVelocityNormal = aTrack.GetVelocity() * 
                                             (OldMomentum * theGlobalNormal);

  G4double Enormal = theMomentumNormal*theMomentumNormal/2./neutron_mass_c2;
  G4double Energy  = aTrack.GetKineticEnergy();
   
  G4double FermiPot2  = 0.;
  G4double pDiffuse   = 0.;
  G4double pSpinFlip  = 0.;
  G4double pUpScatter = 0.;

  if (aMaterialPropertiesTable2) {
     FermiPot2  = aMaterialPropertiesTable2->GetConstProperty("FERMIPOT")*neV;
     pDiffuse   = aMaterialPropertiesTable2->GetConstProperty("DIFFUSION");
     pSpinFlip  = aMaterialPropertiesTable2->GetConstProperty("SPINFLIP");
     pUpScatter = aMaterialPropertiesTable2->GetConstProperty("LOSS");
  }

  G4double FermiPot1 = 0.;
  if (aMaterialPropertiesTable1)
     FermiPot1 = aMaterialPropertiesTable1->GetConstProperty("FERMIPOT")*neV;

  G4double FermiPotDiff = FermiPot2 - FermiPot1;

  if ( verboseLevel > 1 )
     G4cout << "UCNBoundaryProcess: new FermiPot: " << FermiPot2/neV 
            << "neV, old FermiPot:" << FermiPot1/neV << "neV" << G4endl;
  
  // Use microroughness diffuse reflection behavior if activated

  DoMicroRoughnessReflection = UseMicroRoughnessReflection;

  G4double theta_i = 0.;

  if (!aMaterialPropertiesTable2) {

     nNoMPT++;
     theStatus = NoMPT;
     if ( verboseLevel > 0 ) BoundaryProcessVerbose();
     DoMicroRoughnessReflection = false;

  } else {

      if (!aMaterialPropertiesTable2->GetMicroRoughnessTable()) {

         nNoMRT++;
         theStatus = NoMRT;
         if ( verboseLevel > 0 ) BoundaryProcessVerbose();

         DoMicroRoughnessReflection = false;
      }

      // Angle theta_in between surface and momentum direction,
      // Phi_in is defined to be 0

      theta_i = OldMomentum.angle(-theGlobalNormal);

      // Checks the MR-conditions

      if (!aMaterialPropertiesTable2-> 
                    ConditionsValid(Energy, FermiPotDiff, theta_i)) {

         nNoMRCondition++;
         theStatus = NoMRCondition;
         if ( verboseLevel > 0 ) BoundaryProcessVerbose();

          DoMicroRoughnessReflection = false;
      }
  }

  G4double MRpDiffuse = 0.;
  G4double MRpDiffuseTrans = 0.;

  // If microroughness is available and active for material in question

  if (DoMicroRoughnessReflection) {

      // Integral probability for non-specular reflection with microroughness

      MRpDiffuse = aMaterialPropertiesTable2->
                                          GetMRIntProbability(theta_i, Energy);

      // Integral probability for non-specular transmission with microroughness

      MRpDiffuseTrans = aMaterialPropertiesTable2->
                                     GetMRIntTransProbability(theta_i, Energy);

      if ( verboseLevel > 1 ) {
         G4cout << "theta: " << theta_i/degree << "degree" << G4endl;
         G4cout << "Energy: " << Energy/neV << "neV"
                << ", Enormal: " << Enormal/neV << "neV" << G4endl;
         G4cout << "FermiPotDiff: " << FermiPotDiff/neV << "neV " << G4endl;
         G4cout << "Reflection_prob: " << MRpDiffuse 
                << ", Transmission_prob: " << MRpDiffuseTrans << G4endl;
      }
  }

  if (!High(Enormal, FermiPotDiff)) {

     // Below critical velocity

     if (verboseLevel > 0) G4cout << "G4UCNBoundaryProcess -> BELOW critical velocity" << G4endl;

     // Loss on reflection

     if (Loss(pUpScatter, theVelocityNormal, FermiPotDiff)) {

        // kill it.
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        aParticleChange.ProposeLocalEnergyDeposit(Energy);

        nAbsorption++;
        theStatus = Absorption;
        if ( verboseLevel > 0 ) BoundaryProcessVerbose();

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
     }

     // spinflips

     if (SpinFlip(pSpinFlip)) {
        nFlip++;
        theStatus = Flip;
        if ( verboseLevel > 0 ) BoundaryProcessVerbose();

        G4ThreeVector NewPolarization = -1. * aParticle->GetPolarization();
        aParticleChange.ProposePolarization(NewPolarization);
     }

     // Reflect from surface

     G4ThreeVector NewMomentum;

     // If microroughness is available and active - do non-specular reflection

     if (DoMicroRoughnessReflection)
        NewMomentum = MRreflect(MRpDiffuse, OldMomentum, theGlobalNormal,
                                Energy, FermiPotDiff);
     else

        // Else do it with the Lambert model as implemented by Peter Fierlinger

        NewMomentum = Reflect(pDiffuse, OldMomentum, theGlobalNormal);

     aParticleChange.ProposeMomentumDirection(NewMomentum);
    
  } else {

     // Above critical velocity

     if (verboseLevel > 0) G4cout << "G4UCNBoundaryProcess -> ABOVE critical velocity" << G4endl;

     // If it is faster than the criticial velocity,
     // there is a probability to be still reflected.
     // This formula is (only) valid for low loss materials

     // If microroughness available and active, do reflection with it

     G4ThreeVector NewMomentum;

     if (DoMicroRoughnessReflection) {

        G4double Enew;

        NewMomentum = 
               MRreflectHigh(MRpDiffuse, MRpDiffuseTrans, 0., OldMomentum,
                             theGlobalNormal, Energy, FermiPotDiff, Enew);

        if (Enew == 0.) {
          aParticleChange.ProposeTrackStatus(fStopAndKill);
          aParticleChange.ProposeLocalEnergyDeposit(Energy);
          return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        } else {
          aParticleChange.ProposeEnergy(Enew);
          aParticleChange.ProposeMomentumDirection(NewMomentum);
          aParticleChange.ProposeVelocity(std::sqrt(2*Enew/neutron_mass_c2)*c_light);
          aParticleChange.ProposeLocalEnergyDeposit(Energy-Enew);
        }

     } else {

        G4double reflectivity = Reflectivity(FermiPotDiff, Enormal);

        if ( verboseLevel > 1 ) G4cout << "UCNBoundaryProcess: reflectivity "
                                       << reflectivity << G4endl;

        if (G4UniformRand() < reflectivity) { 

          // Reflect from surface

          NewMomentum = Reflect(pDiffuse, OldMomentum, theGlobalNormal);
          aParticleChange.ProposeMomentumDirection(NewMomentum);

        } else {

          // --- Transmission because it is faster than the critical velocity 

          G4double Enew = Transmit(FermiPotDiff, Energy);

          // --- Change of the normal momentum component
          //     p = sqrt(2*m*Ekin)

          G4double mass = -std::sqrt(theMomentumNormal*theMomentumNormal - 
                                neutron_mass_c2*2.*FermiPotDiff);

          // --- Momentum direction in new media

          NewMomentum = 
               theNeutronMomentum - (theMomentumNormal-mass)*theGlobalNormal;

          nSnellTransmit++;
          theStatus = SnellTransmit;
          if ( verboseLevel > 0 ) BoundaryProcessVerbose();

          aParticleChange.ProposeEnergy(Enew);
          aParticleChange.ProposeMomentumDirection(NewMomentum.unit());
          aParticleChange.ProposeVelocity(std::sqrt(2*Enew/neutron_mass_c2)*c_light);
          aParticleChange.ProposeLocalEnergyDeposit(Energy-Enew);

          if (verboseLevel > 1) { 
             G4cout << "Energy: " << Energy/neV << "neV, Enormal: " 
                    << Enormal/neV << "neV, fpdiff " << FermiPotDiff/neV 
                    << "neV, Enew " << Enew/neV << "neV" << G4endl;
	     G4cout << "UCNBoundaryProcess: Transmit and set the new Energy "
                    << aParticleChange.GetEnergy()/neV << "neV" << G4endl;
          }
        }
     }
  }

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4double G4UCNBoundaryProcess::GetMeanFreePath(const G4Track&, 
                                               G4double ,
                                               G4ForceCondition* condition)
{
  *condition = Forced;

  return DBL_MAX;
}

G4bool G4UCNBoundaryProcess::Loss(G4double pUpScatter,
                                  G4double theVelocityNormal,
                                  G4double FermiPot)
{
  // The surface roughness is not taken into account here.
  // One could use e.g. ultracold neutrons, R. Golub, p.35,
  // where mu is increased by roughness parameters sigma and
  // omega, which are related to the height and the distance of
  // "disturbances" on the surface

  G4double vBound = std::sqrt(2.*FermiPot/neutron_mass_c2*c_squared);
  G4double vRatio = theVelocityNormal/vBound;

  G4double pLoss = (2*pUpScatter*vRatio)/(std::sqrt(1-(vRatio*vRatio)));

  // Check, if enhancement for surface roughness should be done

  if (DoMicroRoughnessReflection) {
     if (aMaterialPropertiesTable2) {
        const G4double hdm = hbar_Planck*c_squared/neutron_mass_c2;
        G4double b = aMaterialPropertiesTable2->GetRMS();
        G4double w = aMaterialPropertiesTable2->GetCorrLen();

        // cf. Golub's book p. 35, eq. 2.103

        pLoss *= std::sqrt(1+2*b*b*vBound*vBound/
                    (hdm*hdm+0.85*hdm*vBound*w+2*vBound*vBound*w*w));
     }
  }

  return (G4UniformRand() <= std::fabs(pLoss));
}

G4bool G4UCNBoundaryProcess::SpinFlip(G4double pSpinFlip)
{
  return (G4UniformRand() <= pSpinFlip);
}

G4double G4UCNBoundaryProcess::Reflectivity(G4double FermiPot, G4double Enormal)
{
  G4double r = (std::sqrt(Enormal) - std::sqrt(Enormal - FermiPot)) /
               (std::sqrt(Enormal) + std::sqrt(Enormal - FermiPot));

  return r*r;
} 

G4ThreeVector G4UCNBoundaryProcess::Reflect(G4double pDiffuse, 
                                            G4ThreeVector OldMomentum,
                                            G4ThreeVector Normal)
{
  G4double PdotN = OldMomentum * Normal;

  G4ThreeVector NewMomentum = OldMomentum - (2.*PdotN)*Normal;
  NewMomentum.unit();

  // Reflect diffuse

  if (NewMomentum == OldMomentum || G4UniformRand() < pDiffuse) {

     NewMomentum = LDiffRefl(Normal);

     bLambertianReflection++;
     theStatus = LambertianReflection;
     if ( verboseLevel > 0 ) BoundaryProcessVerbose();

     return NewMomentum;
  }

  // Reflect specular

  bSpecularReflection++;
  theStatus = SpecularReflection;
  if ( verboseLevel > 0 ) BoundaryProcessVerbose();

  return NewMomentum;
}

G4ThreeVector G4UCNBoundaryProcess::MRreflect(G4double pDiffuse,
                                              G4ThreeVector OldMomentum,
                                              G4ThreeVector Normal,
                                              G4double Energy,
                                              G4double FermiPot)
{
  // Only for Enormal <= VFermi

  G4ThreeVector NewMomentum;

  // Checks if the reflection should be non-specular

  if (G4UniformRand() <= pDiffuse) {

      // Reflect diffuse

      // Determines the angles of the non-specular reflection

      NewMomentum =
             MRDiffRefl(Normal, Energy, FermiPot, OldMomentum, pDiffuse);

      bMRDiffuseReflection++;
      theStatus = MRDiffuseReflection;
      if ( verboseLevel > 0 ) BoundaryProcessVerbose();

      return NewMomentum;

  } else {

      // Reflect specular

      G4double PdotN = OldMomentum * Normal;

      NewMomentum = OldMomentum - (2.*PdotN)*Normal;
      NewMomentum.unit();

      bSpecularReflection++;
      theStatus = SpecularReflection;
      if ( verboseLevel > 0 ) BoundaryProcessVerbose();

      return NewMomentum;
  }
}

G4ThreeVector G4UCNBoundaryProcess::MRreflectHigh(G4double pDiffuse,
                                                  G4double pDiffuseTrans,
                                                  G4double pLoss,
                                                  G4ThreeVector OldMomentum,
                                                  G4ThreeVector Normal,
                                                  G4double Energy,
                                                  G4double FermiPot,
                                                  G4double &Enew)
{
  // Only for Enormal > VFermi

  G4double costheta = OldMomentum*Normal;

  G4double Enormal = Energy * (costheta*costheta);

//  G4double pSpecular = Reflectivity(Enormal,FermiPot)*
  G4double pSpecular = Reflectivity(FermiPot,Enormal)*
                                         (1.-pDiffuse-pDiffuseTrans-pLoss);

  G4ThreeVector NewMomentum;

  G4double decide = G4UniformRand();

  if (decide < pSpecular) {

     // Reflect specularly

     G4double PdotN = OldMomentum * Normal;
     NewMomentum = OldMomentum - (2.*PdotN)*Normal;
     NewMomentum.unit();

     Enew = Energy;

     aSpecularReflection++;
     theStatus = SpecularReflection;
     if ( verboseLevel ) BoundaryProcessVerbose();

     return NewMomentum;
  }

  if (decide < pSpecular + pDiffuse) {

     // Reflect diffusely

      // Determines the angles of the non-specular reflection

      NewMomentum =
          MRDiffRefl(Normal, Energy, FermiPot, OldMomentum, pDiffuse);

      if (verboseLevel > 0) G4cout << "Diffuse normal " << Normal
                                   << ", " << NewMomentum << G4endl;
      Enew = Energy;

      aMRDiffuseReflection++;
      theStatus = MRDiffuseReflection;
      if ( verboseLevel ) BoundaryProcessVerbose();

      return NewMomentum;
  }

  if (decide < pSpecular + pDiffuse + pDiffuseTrans) {

     // Transmit diffusely

     // Determines the angles of the non-specular transmission

     NewMomentum =
      MRDiffTrans(Normal, Energy, FermiPot, OldMomentum, pDiffuseTrans);

     Enew = Energy - FermiPot;

     aMRDiffuseTransmit++;
     theStatus = MRDiffuseTransmit;
     if ( verboseLevel ) BoundaryProcessVerbose();

     return NewMomentum;
  }

  if (decide < pSpecular + pDiffuse + pDiffuseTrans + pLoss) {

     // Loss

     Enew = 0.;

     nEzero++;
     theStatus = Ezero;
     if ( verboseLevel > 0 ) BoundaryProcessVerbose();

     return NewMomentum;
  }

  // last case: Refractive transmission

  Enew = Energy - FermiPot;

  G4double palt = std::sqrt(2*neutron_mass_c2/c_squared*Energy);
  G4double produ = OldMomentum * Normal;

  NewMomentum = palt*OldMomentum-
                (palt*produ+std::sqrt(palt*palt*produ*produ-2*neutron_mass_c2/
                      c_squared*FermiPot))*Normal;

  mSnellTransmit++;
  theStatus = SnellTransmit;
  if ( verboseLevel > 0 ) BoundaryProcessVerbose();

  return NewMomentum.unit();
}

G4ThreeVector G4UCNBoundaryProcess::MRDiffRefl(G4ThreeVector Normal,
                                               G4double Energy,
                                               G4double FermiPot,
                                               G4ThreeVector OldMomentum,
                                               G4double pDiffuse)
{
  G4bool accepted = false;

  G4double theta_o, phi_o;

  // Polar angle of incidence

  G4double theta_i = OldMomentum.polarAngle(-Normal);

  // Azimuthal angle of incidence

  //  G4double phi_i = -OldMomentum.azimAngle(-Normal);

  // accept-reject method for MR-reflection

  G4int count = 0;
  while (!accepted) {
        theta_o = G4UniformRand()*pi/2;
        phi_o = G4UniformRand()*pi*2-pi;
        // Box over distribution is increased by 50% to ensure no value is above
        if (1.5*G4UniformRand()*
           aMaterialPropertiesTable2->
             GetMRMaxProbability(theta_i, Energy)/pDiffuse <=
           aMaterialPropertiesTable2->
             GetMRProbability(theta_i,Energy,FermiPot,theta_o,phi_o)/pDiffuse)

           accepted = true;

         // For the case that the box is nevertheless exceeded

        if (aMaterialPropertiesTable2->
             GetMRProbability(theta_i, Energy, FermiPot, theta_o, phi_o)/
             (1.5*aMaterialPropertiesTable2->
                              GetMRMaxProbability(theta_i, Energy)) > 1) {
            G4cout << "MRMax Wahrscheinlichkeitsueberschreitung!" << G4endl;
            G4cout << aMaterialPropertiesTable2->
                   GetMRProbability(theta_i, Energy, FermiPot, theta_o, phi_o)/
                   (1.5*aMaterialPropertiesTable2->
                               GetMRMaxProbability(theta_i, Energy)) << G4endl;
            aMaterialPropertiesTable2->
               SetMRMaxProbability(theta_i, Energy,
                                   aMaterialPropertiesTable2->
                                    GetMRProbability(theta_i, Energy, 
                                                     FermiPot, theta_o, phi_o));
        }
	// Loop checking, 31-Aug-2015, Vladimir Ivanchenko
	if(++count > 10000) { accepted = true; }
  }

  // Creates vector in the local coordinate system of the reflection

  G4ThreeVector localmomentum;
  localmomentum.setRThetaPhi(1., theta_o, phi_o);

  ftheta_o = theta_o;
  fphi_o = phi_o;

  // Get coordinate transform matrix

  G4RotationMatrix TransCoord = 
      GetCoordinateTransformMatrix(Normal, OldMomentum);

  //transfer to the coordinates of the global system

  G4ThreeVector momentum = TransCoord*localmomentum;

  //momentum.rotateUz(Normal);

  if (momentum * Normal<0) {
     momentum*=-1;
     // something has gone wrong...
     G4cout << "G4UCNBoundaryProcess::MRDiffRefl: !" << G4endl;
  }

  return momentum.unit();
}

G4ThreeVector G4UCNBoundaryProcess::MRDiffTrans(G4ThreeVector Normal,
                                                G4double Energy,
                                                G4double FermiPot,
                                                G4ThreeVector OldMomentum,
                                                G4double pDiffuseTrans)
{
  G4bool accepted = false;

  G4double theta_o, phi_o;

  // Polar angle of incidence

  G4double theta_i = OldMomentum.polarAngle(-Normal);

  // azimuthal angle of incidence

  //  G4double phi_i = -OldMomentum.azimAngle(-Normal);

  G4int count = 0;
  while (!accepted) {
    theta_o = G4UniformRand()*pi/2;
    phi_o = G4UniformRand()*pi*2-pi;

    // Box over distribution is increased by 50% to ensure no value is above

    if (1.5*G4UniformRand()*
        aMaterialPropertiesTable2->
          GetMRMaxTransProbability(theta_i, Energy)/pDiffuseTrans <=
        aMaterialPropertiesTable2->
          GetMRTransProbability(theta_i,Energy,FermiPot,theta_o,phi_o)/
                                                          pDiffuseTrans)

        accepted=true;

    // For the case that the box is nevertheless exceeded

    if(aMaterialPropertiesTable2->
        GetMRTransProbability(theta_i, Energy, FermiPot, theta_o, phi_o)/
        (1.5*aMaterialPropertiesTable2->
                         GetMRMaxTransProbability(theta_i, Energy)) > 1) {
        G4cout << "MRMaxTrans Wahrscheinlichkeitsueberschreitung!" << G4endl;
        G4cout << aMaterialPropertiesTable2->
               GetMRTransProbability(theta_i, Energy, FermiPot, theta_o, phi_o)/
               (1.5*aMaterialPropertiesTable2->
                           GetMRMaxTransProbability(theta_i, Energy)) << G4endl;
        aMaterialPropertiesTable2->
           SetMRMaxTransProbability(theta_i, Energy,
                               aMaterialPropertiesTable2->
                                GetMRTransProbability(theta_i, Energy,
                                                 FermiPot, theta_o, phi_o));
    }
    // Loop checking, 31-Aug-2015, Vladimir Ivanchenko
    if(++count > 10000) { accepted = true; }
  }

  // Creates vector in the local coordinate system of the reflection

  G4ThreeVector localmomentum;
  localmomentum.setRThetaPhi(1., pi-theta_o, phi_o);

  // Get coordinate transform matrix

  G4RotationMatrix TransCoord = 
    GetCoordinateTransformMatrix(Normal, OldMomentum);

  // Transfer to the coordinates of the global system

  G4ThreeVector momentum = TransCoord*localmomentum;

  if (momentum*Normal<0) {
     // something has gone wrong... 
     momentum*=-1;
     G4cout << "G4UCNBoundaryProcess::MRDiffTrans: !" << G4endl;
  }

  return momentum.unit();
}

G4double G4UCNBoundaryProcess::Transmit(G4double FermiPot, G4double Energy)
{
  return (Energy - FermiPot);
}

G4ThreeVector G4UCNBoundaryProcess::LDiffRefl(G4ThreeVector Normal)
{
  G4ThreeVector momentum;

  // cosine distribution - Lambert's law

  momentum.setRThetaPhi(1., std::acos(std::sqrt(G4UniformRand())), 2.*pi*G4UniformRand());
  momentum.rotateUz(Normal);

  if (momentum*Normal < 0) {
     momentum*=-1;
     G4cout << "G4UCNBoundaryProcess::LDiffRefl: !" << G4endl;
  }

  return momentum.unit();
}

G4RotationMatrix G4UCNBoundaryProcess::
                GetCoordinateTransformMatrix(G4ThreeVector Normal,
                                             G4ThreeVector direction)
{
   // Definition of the local coordinate system (c.s. of the reflection)

  G4ThreeVector  es1, es2, es3;

  // z-Coordinate is the surface normal, has already length 1

  es3 = Normal;

  // Perpendicular part of incident direction w.r.t. normal
  es1 = direction.perpPart(Normal);

  // Set to unit length: x-Coordinate

  es1.setMag(1.);
  es2 = es1;

  // y-Coordinate is the pi/2-rotation of x-axis around z-axis

  es2.rotate(90.*degree, es3);

  // Transformation matrix consists just of the three coordinate vectors
  // as the global coordinate system is the 'standard' coordinate system

  G4RotationMatrix matrix(es1, es2, es3);

  return matrix;
}

void G4UCNBoundaryProcess::BoundaryProcessVerbose() const
{
  if ( theStatus == Undefined )
     G4cout << " *** Undefined *** " << G4endl;
  if ( theStatus == NotAtBoundary )
     G4cout << " *** NotAtBoundary *** " << G4endl;
  if ( theStatus == SameMaterial )
     G4cout << " *** SameMaterial *** " << G4endl;
  if ( theStatus == StepTooSmall )
     G4cout << " *** StepTooSmall *** " << G4endl;
  if ( theStatus == NoMPT )
     G4cout << " *** No G4UCNMaterialPropertiesTable *** " << G4endl;
  if ( theStatus == NoMRT )
     G4cout << " *** No MicroRoughness Table *** " << G4endl;
  if ( theStatus == NoMRCondition )
     G4cout << " *** MicroRoughness Condition not satisfied *** " << G4endl;
  if ( theStatus == Absorption )
     G4cout << " *** Loss on Surface *** " << G4endl;
  if ( theStatus == Ezero )
     G4cout << " *** Ezero on Surface *** " << G4endl;
  if ( theStatus == Flip )
     G4cout << " *** Spin Flip on Surface *** " << G4endl;
  if ( theStatus == SpecularReflection )
     G4cout << " *** Specular Reflection *** " << G4endl;
  if ( theStatus == LambertianReflection )
     G4cout << " *** LambertianR Reflection *** " << G4endl;
  if ( theStatus == MRDiffuseReflection )
     G4cout << " *** MR Model Diffuse Reflection *** " << G4endl;
  if ( theStatus == SnellTransmit )
     G4cout << " *** Snell Transmission *** " << G4endl;
  if ( theStatus == MRDiffuseTransmit )
     G4cout << " *** MR Model Diffuse Transmission *** " << G4endl;
}

void G4UCNBoundaryProcess::BoundaryProcessSummary(void) const
{
  G4cout << "Sum NoMT:                            "
         << nNoMPT << G4endl;
  G4cout << "Sum NoMRT:                           "
         << nNoMRT << G4endl;
  G4cout << "Sum NoMRCondition:                   "
         << nNoMRCondition << G4endl;
  G4cout << "Sum No. E < V Loss:                  "
         << nAbsorption << G4endl;
  G4cout << "Sum No. E > V Ezero:                 "
         << nEzero << G4endl;
  G4cout << "Sum No. E < V SpinFlip:              "
         << nFlip << G4endl;
  G4cout << "Sum No. E > V Specular Reflection:   "
         << aSpecularReflection << G4endl;
  G4cout << "Sum No. E < V Specular Reflection:   "
         << bSpecularReflection << G4endl;
  G4cout << "Sum No. E < V Lambertian Reflection: "
         << bLambertianReflection << G4endl;
  G4cout << "Sum No. E > V MR DiffuseReflection:  "
         << aMRDiffuseReflection << G4endl;
  G4cout << "Sum No. E < V MR DiffuseReflection:  "
         << bMRDiffuseReflection << G4endl;
  G4cout << "Sum No. E > V SnellTransmit:         "
         << nSnellTransmit << G4endl;
  G4cout << "Sum No. E > V MR SnellTransmit:      "
         << mSnellTransmit << G4endl;
  G4cout << "Sum No. E > V DiffuseTransmit:       "
         << aMRDiffuseTransmit << G4endl;
  G4cout << "                                     " << G4endl;
}

G4bool G4UCNBoundaryProcess::InvokeSD(const G4Step* pStep)
{
  G4Step aStep = *pStep;

  aStep.AddTotalEnergyDeposit(pStep->GetTrack()->GetKineticEnergy());

  G4VSensitiveDetector* sd = aStep.GetPostStepPoint()->GetSensitiveDetector();
  if (sd) return sd->Hit(&aStep);
  else return false;
}
