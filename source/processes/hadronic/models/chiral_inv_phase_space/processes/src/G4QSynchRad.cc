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
// $Id$
//
// Created by Mikhail Kosov 6-Nov-2009
//
// --------------------------------------------------------------
// Short description: Algorithm of Synchrotron Radiation from PDG
// gamma>>1: dI/dw=(8pi/9)*alpha*gamma*F(w/wc), wc=3*gamma^3*c/2/R
// F(y)=(9*sqrt(3)/8/pi)*y*int{y,inf}(K_(5/3)(x)dx) (approximated)
// N_gamma=[5pi/sqrt(3)]*alpha*gamma; <w>=[8/15/sqrt(3)]*wc
// for electrons/positrons: wc(keV)=2.22*[E(GeV)]^3/R(m)
// dE per revolution = (4pi/3)*e^2*beta^3*gamma/R
// at beta=1, dE(MeV)=.o885*[E(GeV)]^4/R(m)
//---------------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QSynchRad.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicDeprecate.hh"


// Constructor
G4QSynchRad::G4QSynchRad(const G4String& Name, G4ProcessType Type):
  G4VDiscreteProcess (Name, Type), minGamma(227.), Polarization(0.,0.,1.) {
  G4HadronicDeprecate("G4QSynchRad");
}

// Calculates MeanFreePath in GEANT4 internal units
G4double G4QSynchRad::GetMeanFreePath(const G4Track& track,G4double,G4ForceCondition* cond)
{
  static const G4double coef = 0.4*std::sqrt(3.)/fine_structure_const;
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  *cond = NotForced ;
  G4double gamma = particle->GetTotalEnergy() / particle->GetMass();
#ifdef debug
  G4cout<<"G4QSynchRad::MeanFreePath: gamma = "<<gamma<<G4endl;
#endif
  G4double MFP = DBL_MAX;
  if( gamma > minGamma )                            // For smalle gamma neglect the process
  {
    G4double R = GetRadius(track);
#ifdef debug
    G4cout<<"G4QSynchRad::MeanFreePath: Radius = "<<R/meter<<" [m]"<<G4endl;
#endif
    if(R > 0.) MFP= coef*R/gamma;
  }
#ifdef debug
  G4cout<<"G4QSynchRad::MeanFreePath = "<<MFP/centimeter<<" [cm]"<<G4endl;
#endif
  return MFP; 
} 

G4VParticleChange* G4QSynchRad::PostStepDoIt(const G4Track& track, const G4Step& step)

{
  static const G4double hc = 1.5 * c_light * hbar_Planck; // E_c=h*w_c=1.5*(hc)*(gamma^3)/R
  aParticleChange.Initialize(track);
  const G4DynamicParticle* particle=track.GetDynamicParticle();
  G4double gamma = particle->GetTotalEnergy() / particle->GetMass();
  if(gamma <= minGamma )
  {
#ifdef debug
    G4cout<<"-Warning-G4QSynchRad::PostStepDoIt is called for small gamma="<<gamma<<G4endl;
#endif
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
  // Photon energy calculation (E < 8.1*Ec restriction)
  G4double R = GetRadius(track);
  if(R <= 0.)
  {
#ifdef debug
    G4cout<<"-Warning-G4QSynchRad::PostStepDoIt: zero or negativ radius ="
          <<R/meter<<" [m]"<<G4endl;
#endif
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }
  G4double EPhoton = hc * gamma * gamma * gamma / R;           // E_c
  G4double dd=5.e-8;
  G4double rnd=G4UniformRand()*(1.+dd);
  if     (rnd < 0.5 ) EPhoton *= .65 * rnd * rnd * rnd;
  else if(rnd > .997) EPhoton *= 15.-1.03*std::log((1.-rnd)/dd+1.);
  else
  {
    G4double r2=rnd*rnd;
    G4double dr=1.-rnd;
    EPhoton*=(2806.+28./rnd)/(1.+500./r2/r2+6500.*(std::sqrt(dr)+28.*dr*dr*dr));
  }
#ifdef debug
  G4cout<<"G4SynchRad::PostStepDoIt: PhotonEnergy = "<<EPhoton/keV<<" [keV]"<<G4endl;
#endif
  if(EPhoton <= 0.)
  {
    G4cout<<"-Warning-G4QSynchRad::PostStepDoIt: zero or negativ photon energy="
          <<EPhoton/keV<<" [keV]"<<G4endl;
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }
  G4double kinEn = particle->GetKineticEnergy();
  G4double newEn = kinEn - EPhoton ;
  if (newEn > 0.)
  {
    aParticleChange.ProposeEnergy(newEn);
    aParticleChange.ProposeLocalEnergyDeposit (0.); 
  } 
  else                                                // Very low probable event
  {
    G4cout<<"-Warning-G4QSynchRad::PostStepDoIt: PhotonEnergy > TotalKinEnergy"<<G4endl;
    EPhoton = kinEn;
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeTrackStatus(fStopButAlive) ;
  } 
  G4ThreeVector MomDir = particle->GetMomentumDirection();
  G4DynamicParticle* Photon = new G4DynamicParticle(G4Gamma::Gamma(), MomDir, EPhoton);
  Photon->SetPolarization(Polarization.x(), Polarization.y(), Polarization.z());
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(Photon);
  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

// Revolution Radius in independent units for the particle (general member function)
G4double G4QSynchRad::GetRadius(const G4Track& track)
{
  static const G4double unk = meter*tesla/0.3/gigaelectronvolt;
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double z = particle->GetDefinition()->GetPDGCharge();
  if(z == 0.) return 0.;                              // --> neutral particle
  if(z < 0.) z=-z;
  G4TransportationManager* transMan = G4TransportationManager::GetTransportationManager();
  G4PropagatorInField* Field = transMan->GetPropagatorInField();
  G4FieldManager* fMan = Field->FindAndSetFieldManager(track.GetVolume());
  if(!fMan || !fMan->GetDetectorField()) return 0.;   // --> no field at all
  const G4Field* pField = fMan->GetDetectorField();
  G4ThreeVector  position = track.GetPosition();
  G4double  PosArray[3]={position.x(), position.y(), position.z()};
  G4double  BArray[3];
  pField->GetFieldValue(PosArray, BArray);
  G4ThreeVector B3D(BArray[0], BArray[1], BArray[2]);
#ifdef debug
  G4cout<<"G4QSynchRad::GetRadius: Pos="<<position/meter<<", B(tesla)="<<B3D/tesla<<G4endl;
#endif
  G4ThreeVector MomDir = particle->GetMomentumDirection();
  G4ThreeVector Ort = B3D.cross(MomDir);
  G4double OrtB = Ort.mag();                          // not negative (independent units)
  if(OrtB == 0.) return 0.;                           // --> along the field line
  Polarization = Ort/OrtB;                            // Polarization unit vector
  G4double mom = particle->GetTotalMomentum();        // Momentum of the particle
#ifdef debug
  G4cout<<"G4QSynchRad::GetRadius: P(GeV)="<<mom/GeV<<", B(tesla)="<<OrtB/tesla<<G4endl;
#endif
  // R [m]= mom [GeV]/(0.3 * z * OrtB [tesla])
  return mom * unk / z / OrtB; 
}
