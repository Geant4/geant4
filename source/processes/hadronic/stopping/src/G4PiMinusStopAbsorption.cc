// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PiMinusStopAbsorption.cc,v 1.3 1999-05-27 09:59:47 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4PiMinusStopAbsorption
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------


#include "G4PiMinusStopAbsorption.hh"
#include <rw/tpordvec.h>
#include <rw/tvordvec.h>
#include <rw/cstring.h>

#include "globals.hh"
#include "Randomize.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTypes.hh"
#include "G4Nucleus.hh"
#include "G4ReactionKinematics.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ThreeVector.hh"

// Constructor

G4PiMinusStopAbsorption::G4PiMinusStopAbsorption(G4PiMinusStopMaterial* materialAlgo, 
						 const G4double Z, const G4double A)
  
{
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
  RWTPtrOrderedVector<G4ParticleDefinition>* defNucleons = _materialAlgo->DefinitionVector();

  G4double newA = _nucleusA;
  G4double newZ = _nucleusZ;

  if (defNucleons != 0)
    {
      for (G4int i=0; i<defNucleons->entries(); i++)
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

  G4double binding = G4NucleiPropertiesTable::GetBindingEnergy(_nucleusZ,_nucleusA) / _nucleusA;
  G4double mass = G4NucleiProperties::GetNuclearMass(newA,newZ);


  RWTPtrOrderedVector<G4LorentzVector>* p4Nucleons = _materialAlgo->P4Vector(binding,mass);

  if (defNucleons != 0 && p4Nucleons != 0)
    {
      int nNucleons = p4Nucleons->entries();
      
      G4double seen = _materialAlgo->FinalNucleons() / 2.;
      G4int maxN = nNucleons;
      if (defNucleons->entries() < nNucleons) { maxN = defNucleons->entries(); }

      for (G4int i=0; i<maxN; i++)
	{
	  G4DynamicParticle* product;
	  if ((*defNucleons)[i] == G4Proton::Proton()) 
	    { product = new G4DynamicParticle(G4Proton::Proton(),*((*p4Nucleons)[i])); }
	  else
	    { product = new G4DynamicParticle(G4Neutron::Neutron(),*((*p4Nucleons)[i])); }
	  G4double ranflat = G4UniformRand();
	  
	  if (ranflat < seen) 
	    { _absorptionProducts->append(product); }
          else
	    { delete product; }
	}
    }

  return _absorptionProducts;

}

G4ThreeVector G4PiMinusStopAbsorption::RecoilMomentum()
{
  G4ThreeVector pProducts(0.,0.,0.);
  
  for (G4int i = 0; i< _absorptionProducts->entries(); i++)
    {
      pProducts = pProducts + (*_absorptionProducts)[i]->GetMomentum();
    }
  return pProducts;
}


G4int G4PiMinusStopAbsorption::NProtons()
{
  G4int n = 0;
  G4int entries = _absorptionProducts->entries();
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
  G4int entries = _absorptionProducts->entries();
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


  G4int nAbsorptionProducts = _absorptionProducts->entries();

  for (int i = 0; i<nAbsorptionProducts; i++)
    {
      productEnergy += (*_absorptionProducts)[i]->GetKineticEnergy();
      pProducts = pProducts + (*_absorptionProducts)[i]->GetMomentum();
      if ((*_absorptionProducts)[i]->GetDefinition() == G4Neutron::Neutron()) nN++;
      if ((*_absorptionProducts)[i]->GetDefinition() == G4Proton::Proton()) nP++;
    }

  G4double productBinding = (G4NucleiPropertiesTable::GetBindingEnergy(_nucleusZ,_nucleusA) / _nucleusA) * nAbsorptionProducts;
  G4double mass = G4NucleiProperties::GetNuclearMass(_nucleusA - (nP + nN),_nucleusZ - nP);
  G4double pNucleus = pProducts.mag();
  G4double eNucleus = sqrt(pNucleus*pNucleus + mass*mass);
  G4double tNucleus = eNucleus - mass;
  G4double temp = G4NucleiPropertiesTable::GetBindingEnergy(_nucleusZ - nP,_nucleusA - (nP + nN)) - 
    G4NucleiPropertiesTable::GetBindingEnergy(_nucleusZ,_nucleusA);
  energy = productEnergy + productBinding + tNucleus;
  
  if (_level > 0)
    {
      cout << "E products " <<  productEnergy  
	   << " Binding " << productBinding << " " << temp << " "
	   << " Tnucleus " << tNucleus 
	   << " energy = " << energy << endl;
    }

  return energy;
}

void G4PiMinusStopAbsorption::SetVerboseLevel(G4int level)
{
  _level = level;
  return;
}
