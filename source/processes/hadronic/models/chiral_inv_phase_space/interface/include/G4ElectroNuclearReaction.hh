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
//
// $Id: G4ElectroNuclearReaction.hh,v 1.15 2003-02-04 10:17:51 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4ElectroNuclearReaction -- header file
// Created: J.P. Wellisch, 12/11/2001
// The last update: J.P. Wellisch, 06-June-02
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
#include "G4GammaParticipants.hh"
#include "G4QGSModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

class G4ElectroNuclearReaction : public G4HadronicInteraction
{
  public: 
    virtual ~G4ElectroNuclearReaction(){}
    G4ElectroNuclearReaction()
    {
      SetMinEnergy(0*GeV);
      SetMaxEnergy(30*TeV);
      
      theHEModel = new G4TheoFSGenerator;
      theCascade = new G4GeneratorPrecompoundInterface;
    }
    
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theLEModel;
    G4TheoFSGenerator * theHEModel;
    G4GeneratorPrecompoundInterface * theCascade;
    G4QGSModel< G4GammaParticipants > theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
    G4ElectroNuclearCrossSection theElectronData;
    G4PhotoNuclearCrossSection thePhotonData;
    G4ParticleChange theResult;
};

inline G4VParticleChange* G4ElectroNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  static const G4double dM=G4Proton::Proton()->GetPDGMass()+G4Neutron::Neutron()->GetPDGMass(); // Mean double nucleon mass = m_n+m_p (@@ no binding)
  static const G4double me=G4Electron::Electron()->GetPDGMass();      // electron mass
  static const G4double me2=me*me;        // squared electron mass
  static const G4double dpi=2*M_PI; // 2*pi
  const G4DynamicParticle* theElectron=aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aD = theElectron->GetDefinition();
  if((aD != G4Electron::ElectronDefinition()) && (aD != G4Positron::PositronDefinition()))
    G4Exception("G4ElectroNuclearReaction::ApplyYourself called for neither electron or positron");
  
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

  // Note: high energy gamma nuclear now implemented.
  
  G4double xSec = theElectronData.GetCrossSection(theElectron, anElement); // Check cross section
  if(xSec<=0.) 
  {
    theResult.Initialize(aTrack);
    return &theResult;        // DO-NOTHING condition
  }
  G4double photonEnergy = theElectronData.GetEquivalentPhotonEnergy();
  G4double theElectronKinEnergy=theElectron->GetKineticEnergy();
  if( theElectronKinEnergy < photonEnergy )
  {
    G4cout << "G4ElectroNuclearReaction::ApplyYourself: photonEnergy is very high"<<G4endl;
    G4cout << "If this condition appears frequently, please contact Hans-Peter.Wellisch@cern.ch"<<G4endl;
  }
  G4double photonQ2 = theElectronData.GetEquivalentPhotonQ2(photonEnergy);
  G4double W=photonEnergy-photonQ2/dM;   // Hadronic energy flow (W-energy) from the virtual photon
  if(getenv("debug_G4ElectroNuclearReaction") )
  {
    G4cout << "G4ElectroNuclearReaction: Equivalent Energy = "<<W<<G4endl;
  }
  if(W<0.) 
  {
    theResult.Initialize(aTrack);
    return &theResult;        // DO-NOTHING condition
    G4Exception("G4ElectroNuclearReaction::ApplyYourself: negative equivalent energy");
  }
  G4DynamicParticle* theDynamicPhoton = new 
                     G4DynamicParticle(G4Gamma::GammaDefinition(), 
		     G4ParticleMomentum(1.,0.,0.), photonEnergy*MeV);            //->-*
  G4double sigNu=thePhotonData.GetCrossSection(theDynamicPhoton, anElement); //       |
  theDynamicPhoton->SetKineticEnergy(W); // Redefine photon with equivalent energy    |
  G4double sigK =thePhotonData.GetCrossSection(theDynamicPhoton, anElement); //       |
  delete theDynamicPhoton; // <-------------------------------------------------------*
  G4double rndFraction = theElectronData.GetVirtualFactor(photonEnergy, photonQ2);
  if(sigNu*G4UniformRand()>sigK*rndFraction) 
  {
    theResult.Initialize(aTrack);
    return &theResult; // DO-NOTHING condition
  }
  theResult.Initialize(aTrack);
  theResult.Clear();
  theResult.SetStatusChange(fAlive);
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
  localTrack.SetStep(aTrack.GetStep());
  G4VParticleChange * result;
  if(photonEnergy < 3*GeV)  
  {
    result = theLEModel.ApplyYourself(localTrack, aTargetNucleus, &theResult);
  } 
  else 
  {
    // G4cout << "0) Getting a high energy electro-nuclear reaction"<<G4endl;
    theHEModel->SetTransport(theCascade);
    theHEModel->SetHighEnergyGenerator(&theStringModel);
    theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
    theStringModel.SetFragmentationModel(theStringDecay);
    theHEModel->SetMinEnergy(2.5*GeV);
    theHEModel->SetMaxEnergy(100*TeV);

    G4VParticleChange * aResult = theHEModel->ApplyYourself(localTrack, aTargetNucleus);
    theResult.SetNumberOfSecondaries(aResult->GetNumberOfSecondaries());
    for(G4int all = 0; all < aResult->GetNumberOfSecondaries(); all++)
    {
      theResult.AddSecondary(aResult->GetSecondary(all));
    }
    result = &theResult;
  }
  return result;
}

#endif
