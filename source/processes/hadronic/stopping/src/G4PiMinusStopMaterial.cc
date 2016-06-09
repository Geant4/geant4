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
//      File name:     G4PiMinusStopMaterial
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
// -------------------------------------------------------------------

#include <vector>

#include "G4PiMinusStopMaterial.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4ParticleTypes.hh"
#include "G4ReactionKinematics.hh"
#include "G4DynamicParticleVector.hh"
#include "G4LorentzVector.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4DistributionGenerator.hh"


// Constructor

G4PiMinusStopMaterial::G4PiMinusStopMaterial()
{
  _definitions = 0;
  _momenta = 0;
  _distributionE = 0;
  _distributionAngle = 0;
  theR = 0.5;
}


// Destructor

G4PiMinusStopMaterial::~G4PiMinusStopMaterial()
{
  if (_definitions != 0) delete _definitions;
  _definitions = 0;

  //A.R. 26-Jul-2012 Coverity fix
  if (_momenta != 0) {
    for (unsigned int i=0; i<_momenta->size(); i++) delete(*_momenta)[i];
    delete _momenta;
  }

  delete _distributionE;
  delete _distributionAngle;
}

std::vector<G4ParticleDefinition*>* G4PiMinusStopMaterial::DefinitionVector()
{

  _definitions->push_back(G4Neutron::Neutron());

  G4double ranflat = G4UniformRand();
  if (ranflat < theR)
    { _definitions->push_back(G4Proton::Proton()); }
  else
    { _definitions->push_back(G4Neutron::Neutron()); }
  
  return _definitions;

}

std::vector<G4LorentzVector*>*
G4PiMinusStopMaterial::P4Vector(const G4double binding,
                                const G4double massNucleus)
{
  // Generate energy of direct absorption products according to experimental
  // data.  The energy distribution of the two nucleons is assumed to be the
  // same for protons and neutrons.  

  G4double eKin1;
  G4double eKin2;
  G4double eRecoil;

  // Assume absorption on two nucleons
  G4int nNucleons = 2;
  G4double availableE = G4PionMinus::PionMinus()->GetPDGMass() - nNucleons * binding;
  G4LorentzVector p1;
  G4LorentzVector p2;

  do 
    {  
      G4double ranflat;
      G4double p;
      G4double energy;
      G4double mass;

      ranflat = G4UniformRand();
      eKin1 = _distributionE->Generate(ranflat);
      mass = (*_definitions)[0]->GetPDGMass();
      energy = eKin1 + mass;
      p = std::sqrt(energy*energy - mass*mass);
      G4double theta1 = pi*G4UniformRand();
      G4double phi1 = GenerateAngle(2.*pi);
      p1 = MakeP4(p,theta1,phi1,energy);

      ranflat = G4UniformRand();
      eKin2 = _distributionE->Generate(ranflat);
      mass = (*_definitions)[1]->GetPDGMass();
      energy = eKin2 + mass;
      p = std::sqrt(energy*energy - mass*mass);
      ranflat = G4UniformRand();
      G4double opAngle = _distributionAngle->Generate(ranflat);
      G4double theta2 = theta1 + opAngle;
      G4double phi2 = phi1 + opAngle;
  
      p2 = MakeP4(p,theta2,phi2,energy);

      G4double pNucleus = (p1.vect() + p2.vect()).mag();
      eRecoil = std::sqrt(pNucleus*pNucleus + massNucleus*massNucleus) - massNucleus;

      // ---- Debug     
      //      G4cout << " ---- binding = " << binding << ", nucleus mass = " << massNucleus 
      //	     << ", p nucleus = " << pNucleus << G4endl;
      //      G4cout << "eKin1,2 " << eKin1 << " " << eKin2 << " eRecoil " << eRecoil 
      //	     << " availableE " << availableE << G4endl;
      // ----

    }  while ((eKin1 + eKin2 + eRecoil) > availableE);
  
  //A.R. 26-Jul-2012 Coverity fix
  if (_momenta != 0) {
    _momenta->push_back(new G4LorentzVector(p1));
    _momenta->push_back(new G4LorentzVector(p2));
  }

  return _momenta;

}

G4double G4PiMinusStopMaterial::GenerateAngle(G4double x)
{
  G4double ranflat = G4UniformRand();
  G4double value = ranflat * x;
  return value;
}

G4LorentzVector G4PiMinusStopMaterial::MakeP4(G4double p, G4double theta, G4double phi, G4double e)
{
  //  G4LorentzVector p4;
  G4double px = p * std::sin(theta) * std::cos(phi);
  G4double py = p * std::sin(theta) * std::sin(phi);
  G4double pz = p * std::cos(theta);
  G4LorentzVector p4(px,py,pz,e);
  return p4;
}

G4double G4PiMinusStopMaterial::RecoilEnergy(const G4double mass)
{
  G4ThreeVector p(0.,0.,0.);

  //A.R. 26-Jul-2012 Coverity fix
  if (_momenta != 0) {  
    for (unsigned int i = 0; i< _momenta->size(); i++)
      {
        p = p + (*_momenta)[i]->vect();
      }
  }

  G4double pNucleus = p.mag();
  G4double eNucleus = std::sqrt(pNucleus*pNucleus + mass*mass);

  return eNucleus;
}

