#include "GhadPotential.hh"
#include "G4NucleiProperties.hh"
#include "G4Nucleon.hh"

// methods for calculating potentials for different types of particles
// aPosition is relative to the nucleus center
G4double GhadPotential::GetPotential(G4Neutron * aN, G4double radius)
{
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);   
   G4double fermiMomentum = theFermi.GetFermiMomentum(density);
   
   const G4double Mn  = aN->GetPDGMass(); 
   G4double result = sqr(fermiMomentum)/(2 * Mn) 
                   + G4NucleiProperties::GetBindingEnergy(theA, theZ)/theA;
   return result;
}

G4double GhadPotential::GetPotential(G4Proton * aP, G4double radius)
{
   G4double theCoulombBarrier = coulomb * theZ/(1. + pow(theA, 1./3.));
   const G4double Mp  = aP->GetPDGMass(); 
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   G4double fermiMomentum = theFermi.GetFermiMomentum(density);

   G4double result = sqr(fermiMomentum)/ (2 * Mp) 
      + G4NucleiProperties::GetBindingEnergy(theA, theZ)/theA 
      + theCoulombBarrier;
   return result;
}

G4double GhadPotential::GetPotential(G4AntiProton * aP, G4double radius)
{
   /*
   G4double theM = theZ * G4Proton::Proton()->GetPDGMass() 
      + (theA - theZ) * G4Neutron::Neutron()->GetPDGMass()
      + G4CreateNucleus::GetBindingEnergy(theZ, theA);
      
   const G4double Mp  = 938.27231 * MeV; // mass of proton
   G4double mu = (theM * Mp)/(theM + Mp);
   
   // antiproton's potential coefficient
   //   V = coeff_antiproton * nucleus_density
   G4double coeff_antiproton = -2.*pi/mu * (1. + Mp) * a_antiproton;
   
   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);
   
   // GetDencity() accepts only G4ThreeVector so build it:   
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;
   
   return coeff_antiproton * density;
   */

   return 0.0;
}

G4double GhadPotential::GetPotential(G4KaonPlus * aK, G4double radius)
{
  return GetKaonPotential(radius);
}
G4double GhadPotential::GetPotential(G4KaonMinus * aK, G4double radius)
{
  return GetKaonPotential(radius);
}
G4double GhadPotential::GetPotential(G4KaonZero * aK, G4double radius)
{
  return GetKaonPotential(radius);
}
G4double GhadPotential::GetPotential(G4AntiKaonZero * aK, G4double radius)
{
  return GetKaonPotential(radius);
}

G4double GhadPotential::GetKaonPotential(G4double radius)
{
   /*
   //G4double theM = G4NucleiProperties::GetAtomicMass(theA, theZ);
   G4double theM = theZ * G4Proton::Proton()->GetPDGMass() 
      + (theA - theZ) * G4Neutron::Neutron()->GetPDGMass()
      + G4CreateNucleus::GetBindingEnergy(theZ, theA);
      
   const G4double Mk  = 496. * MeV;      // mass of "kaon"
   G4double mu = (theM * Mk)/(theM + Mk);
   
   // kaon's potential coefficient
   //   V = coeff_kaon * nucleus_density
   G4double coeff_kaon = -2.*pi/mu * (1. + Mk/theM) * a_kaon;
   
   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);
   
   // GetDencity() accepts only G4ThreeVector so build it:   
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;
   
   return coeff_kaon * density;
   */

   return 0.0;
}

G4double GhadPotential::GetPotential(G4PionMinus * aPi, G4double radius)
{
  return GetPionPotential(radius);
}
G4double GhadPotential::GetPotential(G4PionPlus * aPi, G4double radius)
{
  return GetPionPotential(radius);
}
G4double GhadPotential::GetPotential(G4PionZero * aPi, G4double radius)
{
  return GetPionPotential(radius);
}

G4double GhadPotential::GetPionPotential(G4double radius)
{
   /*
   //G4double theM = G4NucleiProperties::GetAtomicMass(theA, theZ);
   G4double theM = theZ * G4Proton::Proton()->GetPDGMass() 
      + (theA - theZ) * G4Neutron::Neutron()->GetPDGMass()
      + G4CreateNucleus::GetBindingEnergy(theZ, theA);
      
   const G4double Mpi = 139. * MeV;      // mass of "pion"
   G4double mu = (theM * Mpi)/(theM + Mpi);

   // pion's potential coefficient
   //   V = coeff_pion * nucleus_density
   G4double coeff_pion = -2.*pi/mu * (1. + Mpi) * a_pion;
   
   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);
   
   // GetDencity() accepts only G4ThreeVector so build it:   
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;
   
   return coeff_pion * density;
   */

   return 0.0;
}
