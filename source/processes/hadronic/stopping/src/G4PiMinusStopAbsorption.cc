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
//      File name:     G4PiMinusStopAbsorption
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
// -------------------------------------------------------------------


#include "G4PiMinusStopAbsorption.hh"
#include <vector>

#include "globals.hh"
#include "Randomize.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTypes.hh"
#include "G4Nucleus.hh"
#include "G4ReactionKinematics.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ThreeVector.hh"
#include "G4HadronicDeprecate.hh"

// Constructor

G4PiMinusStopAbsorption::G4PiMinusStopAbsorption(G4PiMinusStopMaterial* materialAlgo, 
						 const G4double Z, const G4double A)
  
{
  G4HadronicDeprecate("G4PiMinusStopAbsorption");
  _materialAlgo = materialAlgo;
  _nucleusZ = Z;
  _nucleusA = A;
  _level = 0;
  _absorptionProducts = new G4DynamicParticleVector();
}

// Destructor

G4PiMinusStopAbsorption::~G4PiMinusStopAbsorption()
{
  // Memory management of materialAlgo needs better thought (MGP)
  delete _materialAlgo;
  // Who owns it? Memory management is not clear... (MGP)
  //  _absorptionProducts->clearAndDestroy();
  delete _absorptionProducts;
}

G4DynamicParticleVector* G4PiMinusStopAbsorption::DoAbsorption()
{
  std::vector<G4ParticleDefinition*>* defNucleons = _materialAlgo->DefinitionVector();

  G4double newA = _nucleusA;
  G4double newZ = _nucleusZ;

  if (defNucleons != 0)
    {
      for (unsigned int i=0; i<defNucleons->size(); i++)
	{
	  if ( (*defNucleons)[i] == G4Proton::Proton())
	    {
	      newA = newA - 1;
	      newZ = newZ - 1;
	    }
	  if ((*defNucleons)[i] == G4Neutron::Neutron()) 
	    { newA = newA - 1; }
	}
    }

  G4double binding = G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(_nucleusA) ,static_cast<G4int>(_nucleusZ)) / _nucleusA;
  G4double mass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(newA),static_cast<G4int>(newZ));


  std::vector<G4LorentzVector*>* p4Nucleons = _materialAlgo->P4Vector(binding,mass);

  if (defNucleons != 0 && p4Nucleons != 0)
    {
      unsigned int nNucleons = p4Nucleons->size();
      
      G4double seen = _materialAlgo->FinalNucleons() / 2.;
      G4int maxN = nNucleons;
      if (defNucleons->size() < nNucleons) { maxN = defNucleons->size(); }

      for (G4int i=0; i<maxN; i++)
	{
	  G4DynamicParticle* product;
	  if ((*defNucleons)[i] == G4Proton::Proton()) 
	    { product = new G4DynamicParticle(G4Proton::Proton(),*((*p4Nucleons)[i])); }
	  else
	    { product = new G4DynamicParticle(G4Neutron::Neutron(),*((*p4Nucleons)[i])); }
	  G4double ranflat = G4UniformRand();
	  
	  if (ranflat < seen) 
	    { _absorptionProducts->push_back(product); }
          else
	    { delete product; }
	}
    }

  return _absorptionProducts;

}

G4ThreeVector G4PiMinusStopAbsorption::RecoilMomentum()
{
  G4ThreeVector pProducts(0.,0.,0.);
  
  for (unsigned int i = 0; i< _absorptionProducts->size(); i++)
    {
      pProducts = pProducts + (*_absorptionProducts)[i]->GetMomentum();
    }
  return pProducts;
}


G4int G4PiMinusStopAbsorption::NProtons()
{
  G4int n = 0;
  G4int entries = _absorptionProducts->size();
  for (int i = 0; i<entries; i++)
    {
      if ((*_absorptionProducts)[i]->GetDefinition() == G4Proton::Proton())
	{ n = n + 1; }
    }
  return n;
}


G4int G4PiMinusStopAbsorption::NNeutrons()
{
  G4int n = 0;
  G4int entries = _absorptionProducts->size();
  for (int i = 0; i<entries; i++)
    {
      if ((*_absorptionProducts)[i]->GetDefinition() == G4Neutron::Neutron())
	{ n = n + 1; }
    }
  return n;
}


G4double G4PiMinusStopAbsorption::Energy()
{
  G4double energy = 0.;
  G4double productEnergy = 0.;
  G4ThreeVector pProducts(0.,0.,0.);
  G4int nN = 0;
  G4int nP = 0;


  G4int nAbsorptionProducts = _absorptionProducts->size();

  for (int i = 0; i<nAbsorptionProducts; i++)
    {
      productEnergy += (*_absorptionProducts)[i]->GetKineticEnergy();
      pProducts = pProducts + (*_absorptionProducts)[i]->GetMomentum();
      if ((*_absorptionProducts)[i]->GetDefinition() == G4Neutron::Neutron()) nN++;
      if ((*_absorptionProducts)[i]->GetDefinition() == G4Proton::Proton()) nP++;
    }

  G4double productBinding = (G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(_nucleusA),static_cast<G4int>(_nucleusZ)) / _nucleusA) * nAbsorptionProducts;
  G4double mass = G4NucleiProperties::GetNuclearMass(_nucleusA - (nP + nN),_nucleusZ - nP);
  G4double pNucleus = pProducts.mag();
  G4double eNucleus = std::sqrt(pNucleus*pNucleus + mass*mass);
  G4double tNucleus = eNucleus - mass;
  G4double temp = 
    G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(_nucleusA - (nP + nN)),static_cast<G4int>(_nucleusZ - nP)) - 
    G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(_nucleusA),static_cast<G4int>(_nucleusZ));
  energy = productEnergy + productBinding + tNucleus;
  
  if (_level > 0)
    {
      std::cout << "E products " <<  productEnergy  
	   << " Binding " << productBinding << " " << temp << " "
	   << " Tnucleus " << tNucleus 
	   << " energy = " << energy << G4endl;
    }

  return energy;
}

void G4PiMinusStopAbsorption::SetVerboseLevel(G4int level)
{
  _level = level;
  return;
}
