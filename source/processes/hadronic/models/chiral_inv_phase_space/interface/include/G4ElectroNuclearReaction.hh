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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
// $Id: G4ElectroNuclearReaction.hh,v 1.8 2002-05-22 11:43:15 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4ElectroNuclearReaction -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 21-May-02
//
#ifndef G4ElectroNuclearReaction_h
#define G4ElectroNuclearReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4ElectroNuclearCrossSection.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

class G4ElectroNuclearReaction : public G4HadronicInteraction
{
  public: 
    virtual ~G4ElectroNuclearReaction(){}
    
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
    G4ElectroNuclearCrossSection theElectronData;
    G4PhotoNuclearCrossSection thePhotonData;
    G4ParticleChange theResult;
};

inline G4VParticleChange* G4ElectroNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  static const G4double dM=938.27+939.57; // Mean double nucleon mass = m_n+m_p (@@ no binding)
  static const G4double me=.5109989;      // electron mass
  static const G4double me2=me*me;        // squared electron mass
  static const G4double dpi=2*3.14159265; // 2*pi
  const G4DynamicParticle* theElectron=aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aD = theElectron->GetDefinition();
  if((aD != G4Electron::ElectronDefinition()) && (aD != G4Positron::PositronDefinition()))
    G4Exception("G4ElectroNuclearReaction::ApplyYourself called for neither electron or positron");
  
  theResult.Initialize(aTrack);

  const G4ElementTable* aTab = G4Element::GetElementTable();
  G4Element * anElement = 0;
  G4int aZ = static_cast<G4int>(aTargetNucleus.GetZ()+.1);
  for(size_t ii=0; ii<aTab->size(); ii++)
  {
    if ( abs((*aTab)[ii]->GetZ()-aZ) < .1)
    {
      anElement = (*aTab)[ii];
      break;
    }
  }
  if(0==anElement) 
  {
    G4cerr<<"***G4ElectroNuclearReaction::ApplyYourself: element with Z="<<aTargetNucleus.GetZ()<<
	  " is not in the element table"<<G4endl; // @@ how to retrieve A or N for the isotop?
    G4Exception("Anomalous element error.");
  }
  // This reduces the energy deposition at high energy. It is better to have 0-order fragmentation.
  //G4double xSec;
  //G4double photonEnergy = 10*GeV;  
  //while(photonEnergy>3.*GeV)
  //{
  //  xSec = theElectronData.GetCrossSection(theElectron, anElement); // Can be before while
  //  photonEnergy = theElectronData.GetEffectivePhotonEnergy();
  //}
  G4double xSec = theElectronData.GetCrossSection(theElectron, anElement); // Check cross section
  if(xSec<=0.) return &theResult;        // DO-NOTHING condition
  G4double photonEnergy = theElectronData.GetEquivalentPhotonEnergy();
  G4double theElectronKinEnergy=theElectron->GetKineticEnergy();
  if( theElectronKinEnergy < photonEnergy )
    G4Exception("G4ElectroNuclearReaction::ApplyYourself: photonEnergy above electron energy");
  G4double photonQ2 = theElectronData.GetEquivalentPhotonQ2(photonEnergy);
  G4double K=photonEnergy-photonQ2/dM;   // Equivalent energy (W-energy) of the virtual photon
  if(K<0.) G4Exception("G4ElectroNuclearReaction::ApplyYourself: negative equivalent energy");
  G4DynamicParticle* theDynamicPhoton = new G4DynamicParticle(G4Gamma::GammaDefinition(), 
															  G4ParticleMomentum(1.,0.,0.),
															  photonEnergy*MeV); //->-*
  G4double sigNu=thePhotonData.GetCrossSection(theDynamicPhoton, anElement); //       |
  theDynamicPhoton->SetKineticEnergy(K); // Redefine photon with equivalent energy    |
  G4double sigK =thePhotonData.GetCrossSection(theDynamicPhoton, anElement); //       |
  delete theDynamicPhoton; // <-------------------------------------------------------*
  G4double rndFraction = theElectronData.GetVirtualFactor(photonEnergy, photonQ2);
  if(sigNu*G4UniformRand()>sigK*rndFraction) return &theResult; // DO-NOTHING condition
  // Scatter an electron and make gamma+A reaction
  G4double iniE=theElectronKinEnergy+me; // Initial total energy of electron
  G4double finE=iniE-photonEnergy;       // Final total energy of electron
  theResult.SetEnergyChange(finE-me);    // Modifies the KINETIC ENERGY (Why not in the name?)
  G4double EEm=iniE*finE-me2;            // Just an intermediate value to avoid "2*"
  G4double iniP=sqrt(iniE*iniE-me2);     // Initial momentum of the electron
  G4double finP=sqrt(finE*finE-me2);     // Final momentum of the electron
  G4double cost=(EEm+EEm-photonQ2)/iniP/finP; // cos(theta) for the electron scattering
  if(cost>1.) cost=1.;
  if(cost<-1.) cost=-1.;
  G4ThreeVector dir=theElectron->GetMomentumDirection(); // Direction of primary electron
  G4ThreeVector ort=dir.orthogonal();    // Not normed orthogonal vector (!) (to dir)
  G4ThreeVector ortx = ort.unit();       // First unit vector orthogonal to the direction
  G4ThreeVector orty = dir.cross(ortx);  // Second unit vector orthoganal to the direction
  G4double sint=sqrt(1.-cost*cost);      // Perpendicular component
  G4double phi=dpi*G4UniformRand();      // phi of scattered electron
  G4double sinx=sint*sin(phi);           // x-component
  G4double siny=sint*cos(phi);           // y-component
  G4ThreeVector findir=cost*dir+sinx*ortx+siny*orty;
  theResult.SetMomentumDirectionChange(findir); // new direction for the electron
  G4ThreeVector photonMomentum=iniP*dir-finP*findir;
  G4DynamicParticle localGamma(G4Gamma::GammaDefinition(), photonEnergy, photonMomentum);
  //G4DynamicParticle localGamma(G4Gamma::GammaDefinition(), photonDirection, photonEnergy);
  //G4DynamicParticle localGamma(G4Gamma::GammaDefinition(), photonLorentzVector);
  G4ThreeVector position(0,0,0);
  G4Track localTrack(&localGamma, 0., position);
  G4VParticleChange* result = theModel.ApplyYourself(localTrack, aTargetNucleus, &theResult);
  return result;
}

#endif
