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
// $Id: G4PiMinusStopMaterial.cc,v 1.8 2001-10-05 16:09:44 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4PiMinusStopMaterial
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#include "G4ios.hh"

#include "G4PiMinusStopMaterial.hh"

#include "g4std/vector"

#include "globals.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4ParticleTypes.hh"
#include "G4ReactionKinematics.hh"
#include "G4DynamicParticleVector.hh"
#include "G4LorentzVector.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4DistributionGenerator.hh"


// Constructor

G4PiMinusStopMaterial::G4PiMinusStopMaterial()
  
{
  _definitions = 0;
  _momenta = 0;
  _distributionE = 0;
  _distributionAngle = 0;

}


// Destructor

G4PiMinusStopMaterial::~G4PiMinusStopMaterial()
{
  //  _definitions->clear();
  if (_definitions != 0) delete _definitions;
  _definitions = 0;

  for(unsigned int i=0; i<_momenta->size(); i++) delete(*_momenta)[i];
  if (_momenta != 0) delete _momenta;

  delete _distributionE;
  delete _distributionAngle;
}

G4std::vector<G4ParticleDefinition*>* G4PiMinusStopMaterial::DefinitionVector()
{

  _definitions->push_back(G4Neutron::Neutron());

  G4double ranflat = G4UniformRand();
  if (ranflat < theR)
    { _definitions->push_back(G4Proton::Proton()); }
  else
    { _definitions->push_back(G4Neutron::Neutron()); }
  
  return _definitions;

}

G4std::vector<G4LorentzVector*>* G4PiMinusStopMaterial::P4Vector(const G4double binding,
								      const G4double massNucleus)
{

  // Generate energy of direct absorption products according to experimental data
  // The energy distribution of the two nucleons is assumed to be the same 
  // for protons and neutrons  


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
      p = sqrt(energy*energy - mass*mass);
      G4double theta1 = GenerateAngle(pi);
      G4double phi1 = GenerateAngle(2.*pi);
      p1 = MakeP4(p,theta1,phi1,energy);

      ranflat = G4UniformRand();
      eKin2 = _distributionE->Generate(ranflat);
      mass = (*_definitions)[1]->GetPDGMass();
      energy = eKin2 + mass;
      p = sqrt(energy*energy - mass*mass);
      ranflat = G4UniformRand();
      G4double opAngle = _distributionAngle->Generate(ranflat);
      G4double theta2 = theta1 + opAngle;
      G4double phi2 = phi1 + opAngle;
  
      p2 = MakeP4(p,theta2,phi2,energy);

      G4double pNucleus = (p1.vect() + p2.vect()).mag();
      eRecoil = sqrt(pNucleus*pNucleus + massNucleus*massNucleus) - massNucleus;

      // ---- Debug     
      //      G4cout << " ---- binding = " << binding << ", nucleus mass = " << massNucleus 
      //	     << ", p nucleus = " << pNucleus << G4endl;
      //      G4cout << "eKin1,2 " << eKin1 << " " << eKin2 << " eRecoil " << eRecoil 
      //	     << " availableE " << availableE << G4endl;
      // ----

    }  while ((eKin1 + eKin2 + eRecoil) > availableE);
  
  _momenta->push_back(new G4LorentzVector(p1));
  _momenta->push_back(new G4LorentzVector(p2));

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
  G4double px = p * sin(theta) * cos(phi);
  G4double py = p * sin(theta) * sin(phi);
  G4double pz = p * cos(theta);
  G4LorentzVector p4(px,py,pz,e);
  return p4;
}

G4double G4PiMinusStopMaterial::RecoilEnergy(const G4double mass)
{
  G4ThreeVector p(0.,0.,0.);
  
  for (unsigned int i = 0; i< _momenta->size(); i++)
    {
      p = p + (*_momenta)[i]->vect();
    }
  G4double pNucleus = p.mag();
  G4double eNucleus = sqrt(pNucleus*pNucleus + mass*mass);

  return eNucleus;
}










