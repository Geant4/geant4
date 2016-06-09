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
// $Id$
//
// Class description:
// G4ElectroNuclearReaction = eA interface to use CHIPS in the G4Hadronic frame
// Created: J.P. Wellisch, following M. Kossov's algorithm. 12/11/2001
// The last update: J.P. Wellisch, 06-June-02
// 17.02.2009 M.Kossov, now it is recommended to use the G4QCollision process
// 10.11.2010 V.Ivanchenko use cross sections by pointer and not by value
// 07.09.2011 M.Kelsey, follow changes to G4HadFinalState interface

#ifndef G4ElectroNuclearReaction_h
#define G4ElectroNuclearReaction_h 1

#include <iostream>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4ElectroNuclearCrossSection.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitedStringDecay.hh"


class G4ElectroNuclearReaction : public G4HadronicInteraction
{
public: 
  G4ElectroNuclearReaction():G4HadronicInteraction("CHIPSElectroNuclear")
  {
    SetMinEnergy(0*CLHEP::GeV);
    SetMaxEnergy(100*CLHEP::TeV);
      
    theHEModel = new G4TheoFSGenerator;
    theCascade = new G4GeneratorPrecompoundInterface;
    theHEModel->SetTransport(theCascade);
    theHEModel->SetHighEnergyGenerator(&theStringModel);
    theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
    theStringModel.SetFragmentationModel(theStringDecay);
    theHEModel->SetMinEnergy(2.5*CLHEP::GeV);
    theHEModel->SetMaxEnergy(100*CLHEP::TeV);
    theElectronData = new G4ElectroNuclearCrossSection;
    thePhotonData = new G4PhotoNuclearCrossSection;
    Description();
  }

  ~G4ElectroNuclearReaction()
  {
    delete theHEModel;
    delete theCascade;
    delete theStringDecay;
    delete theElectronData; 
    delete thePhotonData;
  }
    
  G4HadFinalState*
  ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);

  void Description() const 
  {
    char* dirName = getenv("G4PhysListDocDir");
    if (dirName) {
      std::ofstream outFile;
      G4String outFileName = GetModelName() + ".html";
      G4String pathName = G4String(dirName) + "/" + outFileName;

      outFile.open(pathName);
      outFile << "<html>\n";
      outFile << "<head>\n";

      outFile << "<title>Description of CHIPS ElectroNuclear Model</title>\n";
      outFile << "</head>\n";
      outFile << "<body>\n";

      outFile << "G4ElectroNuclearReaction handles the inelastic scattering\n"
              << "of e- and e+ from nuclei using the Chiral Invariant Phase\n"
              << "Space (CHIPS) model of M. Kossov.  This model uses the\n"
              << "Equivalent Photon Approximation in which the incoming\n"
              << "electron generates a virtual photon at the electromagnetic\n"
              << "vertex, and the virtual photon is converted to a real photon\n"
              << "before it interacts with the nucleus.  The real photon\n"
              << "interacts with the hadrons in the target by producing\n"
              << "quasmons (or generalized excited hadrons) which then decay\n"
              << "into final state hadrons.  This model is valid for e- and\n"
              << "e+ of all incident energies.\n";

      outFile << "</body>\n";
      outFile << "</html>\n";
      outFile.close();
    }
  }

private:

  G4ChiralInvariantPhaseSpace      theLEModel;
  G4TheoFSGenerator*               theHEModel;
  G4GeneratorPrecompoundInterface* theCascade;
  G4QGSModel<G4GammaParticipants>  theStringModel;
  G4QGSMFragmentation              theFragmentation;
  G4ExcitedStringDecay*            theStringDecay;
  G4ElectroNuclearCrossSection*    theElectronData;
  G4PhotoNuclearCrossSection*      thePhotonData;
  G4HadFinalState                  theResult;
};

inline
G4HadFinalState* G4ElectroNuclearReaction::ApplyYourself(const G4HadProjectile& aTrack, 
                                                               G4Nucleus& aTargetNucleus)
{
  theResult.Clear();
  static const G4double dM=G4Proton::Proton()->GetPDGMass() +
    G4Neutron::Neutron()->GetPDGMass(); // MeanDoubleNucleon Mass = m_n+m_p (@@ no binding)
  static const G4double me=G4Electron::Electron()->GetPDGMass(); // electron mass
  static const G4double me2=me*me;                               // squared electron mass
  G4DynamicParticle theTempEl(const_cast<G4ParticleDefinition *>(aTrack.GetDefinition()), 
                              aTrack.Get4Momentum().vect());
  const G4DynamicParticle* theElectron=&theTempEl;
  const G4ParticleDefinition* aD = theElectron->GetDefinition();
  if((aD != G4Electron::ElectronDefinition()) && (aD != G4Positron::PositronDefinition()))
    throw G4HadronicException(__FILE__, __LINE__,
        "G4ElectroNuclearReaction::ApplyYourself called for neither electron or positron");
  const G4ElementTable* aTab = G4Element::GetElementTable();
  G4Element * anElement = 0;
  G4int aZ = static_cast<G4int>(aTargetNucleus.GetZ_asInt()+.1);
  for(size_t ii=0; ii<aTab->size(); ++ii) if ( std::abs((*aTab)[ii]->GetZ()-aZ) < .1)
  {
    anElement = (*aTab)[ii];
    break;
  }
  if(0==anElement) 
  {
    G4cerr<<"***G4ElectroNuclearReaction::ApplyYourself: element with Z="
          <<aTargetNucleus.GetZ_asInt()<<" is not in the element table"<<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, "Anomalous element error.");
  }

  // Note: high energy gamma nuclear now implemented.
  G4double xSec = theElectronData->GetCrossSection(theElectron, anElement);// Check out XS
  if(xSec<=0.) 
  {
    theResult.SetStatusChange(isAlive);
    theResult.SetEnergyChange(theElectron->GetKineticEnergy());
    // new direction for the electron
    theResult.SetMomentumChange(theElectron->GetMomentumDirection());
    return &theResult;        // DO-NOTHING condition
  }
  G4double photonEnergy = theElectronData->GetEquivalentPhotonEnergy();
  G4double theElectronKinEnergy=theElectron->GetKineticEnergy();
  if( theElectronKinEnergy < photonEnergy )
  {
    G4cout<<"G4ElectroNuclearReaction::ApplyYourself: photonEnergy is very high"<<G4endl;
    G4cout<<">>> If this condition persists, please contact Geant4 group"<<G4endl;
    theResult.SetStatusChange(isAlive);
    theResult.SetEnergyChange(theElectron->GetKineticEnergy());
    // new direction for the electron
    theResult.SetMomentumChange(theElectron->GetMomentumDirection());
    return &theResult;        // DO-NOTHING condition
  }
  G4double photonQ2 = theElectronData->GetEquivalentPhotonQ2(photonEnergy);
  G4double W=photonEnergy-photonQ2/dM; // Hadronic energy flow from the virtual photon
  if(getenv("debug_G4ElectroNuclearReaction") )
  {
    G4cout << "G4ElectroNuclearReaction: Equivalent Energy = "<<W<<G4endl;
  }
  if(W<0.) 
  {
    theResult.SetStatusChange(isAlive);
    theResult.SetEnergyChange(theElectron->GetKineticEnergy());
    // new direction for the electron
    theResult.SetMomentumChange(theElectron->GetMomentumDirection());
    return &theResult;        // DO-NOTHING condition
    // throw G4HadronicException(__FILE__, __LINE__,
    //             "G4ElectroNuclearReaction::ApplyYourself: negative equivalent energy");
  }
  G4DynamicParticle* theDynamicPhoton = new 
                     G4DynamicParticle(G4Gamma::GammaDefinition(), 
                     G4ParticleMomentum(1.,0.,0.), photonEnergy*CLHEP::MeV);         //----->-*
  G4double sigNu=thePhotonData->GetCrossSection(theDynamicPhoton, anElement); //       |
  theDynamicPhoton->SetKineticEnergy(W);  // Redefine photon with equivalent energy    |
  G4double sigK =thePhotonData->GetCrossSection(theDynamicPhoton, anElement); //       |
  delete theDynamicPhoton; // <--------------------------------------------------------*
  G4double rndFraction = theElectronData->GetVirtualFactor(photonEnergy, photonQ2);
  if(sigNu*G4UniformRand()>sigK*rndFraction) 
  {
    theResult.SetStatusChange(isAlive);
    theResult.SetEnergyChange(theElectron->GetKineticEnergy());
    // new direction for the electron
    theResult.SetMomentumChange(theElectron->GetMomentumDirection());
    return &theResult; // DO-NOTHING condition
  }
  theResult.SetStatusChange(isAlive);
  // Scatter an electron and make gamma+A reaction
  G4double iniE=theElectronKinEnergy+me; // Initial total energy of electron
  G4double finE=iniE-photonEnergy;       // Final total energy of electron
  theResult.SetEnergyChange(std::max(0.,finE-me)); // Modifies the KINETIC ENERGY
  G4double EEm=iniE*finE-me2;            // Just an intermediate value to avoid "2*"
  G4double iniP=std::sqrt(iniE*iniE-me2);          // Initial momentum of the electron
  G4double finP=std::sqrt(finE*finE-me2);          // Final momentum of the electron
  G4double cost=(EEm+EEm-photonQ2)/iniP/finP;// std::cos(theta) for the electron scattering
  if(cost> 1.) cost= 1.;
  if(cost<-1.) cost=-1.;
  G4ThreeVector dir=theElectron->GetMomentumDirection(); // Direction of primary electron
  G4ThreeVector ort=dir.orthogonal();    // Not normed orthogonal vector (!) (to dir)
  G4ThreeVector ortx = ort.unit();       // First unit vector orthogonal to the direction
  G4ThreeVector orty = dir.cross(ortx);  // Second unit vector orthoganal to the direction
  G4double sint=std::sqrt(1.-cost*cost);      // Perpendicular component
  G4double phi=CLHEP::twopi*G4UniformRand();      // phi of scattered electron
  G4double sinx=sint*std::sin(phi);           // x-component
  G4double siny=sint*std::cos(phi);           // y-component
  G4ThreeVector findir=cost*dir+sinx*ortx+siny*orty;
  theResult.SetMomentumChange(findir); // new direction for the electron
  G4ThreeVector photonMomentum=iniP*dir-finP*findir;
  G4DynamicParticle localGamma(G4Gamma::GammaDefinition(), photonEnergy, photonMomentum);
  //G4DynamicParticle localGamma(G4Gamma::GammaDefinition(),photonDirection, photonEnergy);
  //G4DynamicParticle localGamma(G4Gamma::GammaDefinition(), photonLorentzVector);
  G4ThreeVector position(0,0,0);
  G4HadProjectile localTrack(localGamma);
  G4HadFinalState * result;
  if (photonEnergy < 3*CLHEP::GeV)
    result = theLEModel.ApplyYourself(localTrack, aTargetNucleus, &theResult);
  else {
    // G4cout << "0) Getting a high energy electro-nuclear reaction"<<G4endl;
    G4HadFinalState * aResult = theHEModel->ApplyYourself(localTrack, aTargetNucleus);
    theResult.AddSecondaries(aResult);
    aResult->Clear();
    result = &theResult;
  }
  return result;
}

#endif
