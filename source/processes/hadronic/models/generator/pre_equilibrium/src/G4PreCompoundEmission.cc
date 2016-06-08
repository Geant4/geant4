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
// $Id: G4PreCompoundEmission.cc,v 1.7 2001/09/27 16:51:56 larazb Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 


#include "G4PreCompoundEmission.hh"

const G4PreCompoundEmission & G4PreCompoundEmission::operator=(const G4PreCompoundEmission &right)
{
  G4Exception("G4PreCompoundEmission::operator= meant to not be accessable");
  return *this;
}


G4bool G4PreCompoundEmission::operator==(const G4PreCompoundEmission &right) const
{
  return false;
}

G4bool G4PreCompoundEmission::operator!=(const G4PreCompoundEmission &right) const
{
  return true;
}


G4PreCompoundEmission::G4PreCompoundEmission(const G4Fragment& aFragment)
{
  // Assume that projectile is a proton
  ProjEnergy = aFragment.GetExcitationEnergy();
  theIncidentDirection = aFragment.GetMomentum().vect().unit();
}


G4ReactionProduct * G4PreCompoundEmission::PerformEmission(G4Fragment & aFragment)
{
  // Choose a Fragment for emission
  G4VPreCompoundFragment * theFragment = theFragmentsVector.ChooseFragment();

  // Kinetic Energy of emitted fragment
  G4double KineticEnergyOfEmittedFragment = theFragment->GetKineticEnergy(aFragment);

  // Calculate the fragment momentum (three vector)
  G4ThreeVector momentum = AngularDistribution(theFragment,aFragment,KineticEnergyOfEmittedFragment);
	
  // Mass of emittef fragment
  G4double EmittedMass = theFragment->GetNuclearMass();

	
  // Now we can calculate the four momentum 
  // both options are valid and give the same result but 2nd one is faster
  //  G4LorentzVector EmittedMomentum(momentum,sqrt(momentum.mag2()+EmittedMass*EmittedMass));
  G4LorentzVector EmittedMomentum(momentum,EmittedMass+KineticEnergyOfEmittedFragment);
	
  // Perform Lorentz boost
  EmittedMomentum.boost(aFragment.GetMomentum().boostVector());  

  // Set emitted fragment momentum
  theFragment->SetMomentum(EmittedMomentum);	


  // NOW THE RESIDUAL NUCLEUS
  // ------------------------

  // Now the residual nucleus. 
  // The energy conservation says that 
  G4double ResidualEcm = 
    aFragment.GetGroundStateMass() + aFragment.GetExcitationEnergy() // initial energy in cm
    - (EmittedMass+KineticEnergyOfEmittedFragment); 

  // Then the four momentum for residual is 
  G4LorentzVector RestMomentum(-momentum,ResidualEcm);
  G4LorentzVector RestMomentum2(aFragment.GetMomentum()-EmittedMomentum);


  // Just for test
  // Excitation energy
  //  G4double anU = ResidualEcm - theFragment->GetRestNuclearMass();
  // This is equivalent
  //  G4double anU = theFragment->GetMaximalKineticEnergy() - KineticEnergyOfEmittedFragment + 
  //    theFragment->GetCoulombBarrier();
	
  // check that Excitation energy is >= 0
  G4double anU = RestMomentum.m()-theFragment->GetRestNuclearMass();
  if (anU < 0.0) G4Exception("G4PreCompoundModel::DeExcite: Excitation energy less than 0!");

  
  
  // Update nucleus parameters:
  // --------------------------
  // Number of excitons
  aFragment.SetNumberOfParticles(aFragment.GetNumberOfParticles()-
				G4int(theFragment->GetA()));
  // Number of charges
  aFragment.SetNumberOfCharged(aFragment.GetNumberOfCharged()-
			       G4int(theFragment->GetZ()));

  // Atomic number
  aFragment.SetA(theFragment->GetRestA());
	  
  // Charge
  aFragment.SetZ(theFragment->GetRestZ());


  // Perform Lorentz boosts
  RestMomentum.boost(aFragment.GetMomentum().boostVector());

  // Update nucleus momentum
  aFragment.SetMomentum(RestMomentum);
	
  // Create a G4ReactionProduct 
  G4ReactionProduct * MyRP = theFragment->GetReactionProduct();
  return MyRP;
}


G4ThreeVector G4PreCompoundEmission::AngularDistribution(G4VPreCompoundFragment * theFragment,
							 const G4Fragment& aFragment,
							 const G4double KineticEnergyOfEmittedFragment) const
{
  G4double p = aFragment.GetNumberOfParticles();
  G4double h = aFragment.GetNumberOfHoles();
  G4double U = aFragment.GetExcitationEnergy();
	
  // Kinetic Energy of emitted fragment
  // G4double KineticEnergyOfEmittedFragment = theFragment->GetKineticEnergy(aFragment);
	
  // Emission particle separation energy
  G4double Bemission = theFragment->GetBindingEnergy();
	
  // Fermi energy
  G4double Ef = G4PreCompoundParameters::GetAddress()->GetFermiEnergy();
	
  //
  G4double g = 0.595*aFragment.GetA()*G4PreCompoundParameters::GetAddress()->GetLevelDensity();
	
  // Average exciton energy relative to bottom of nuclear well
  G4double Eav = 2.0*p*(p+1.0)/((p+h)*g);
	
  // Excitation energy relative to the Fermi Level
  //	G4double Uf = U - (p - h)*Ef;
  G4double Uf = U - KineticEnergyOfEmittedFragment - Bemission;

	
  Eav *= rho(p+1,h,g,Uf,Ef)/rho(p,h,g,Uf,Ef);
	
  Eav += - Uf/(p+h) + Ef;
	
  G4double zeta = G4std::max(1.0,9.3/sqrt(KineticEnergyOfEmittedFragment/MeV));
	
  G4double an = 3.0*sqrt((ProjEnergy+Ef)*(KineticEnergyOfEmittedFragment+Bemission+Ef))/
    (zeta*2.0*aFragment.GetNumberOfExcitons()*Eav);
//  			(zeta*(aFragment.GetNumberOfExcitons()-1.0)*Eav);
			
			
  //  G4double normalization = pi*bessi0(an);
  // Normalization is expesive to calculate
  // And rejection method is selfnormalized
  // So, I only need the maximum
  G4double normalization = exp(an);
	
  G4double theta = 0.0;
  G4double distrib = 0.0;
  do {
    theta = pi*G4UniformRand();
    //    distrib = exp(an*cos(theta))/normalization;
    distrib = exp(an*cos(theta));
  } while ( G4UniformRand()*normalization > distrib );
	
  G4double phi = twopi*G4UniformRand();
	
  // Calculate the momentum magnitude of emitted fragment 	
  G4double EmittedMass = theFragment->GetNuclearMass();
  G4double pmag = sqrt(KineticEnergyOfEmittedFragment*(KineticEnergyOfEmittedFragment+2.0*EmittedMass));
  

  G4double sinTheta = sin(theta);
  //  G4double cosTheta = sqrt(1.0-sinTheta*sinTheta);
  G4double cosTheta = cos(theta);

  G4ThreeVector momentum(pmag*cos(phi)*sinTheta,pmag*sin(phi)*sinTheta,pmag*cosTheta);
  // theta is the angle wrt the incident direction
  momentum.rotateUz(theIncidentDirection);

  return momentum;
}


G4double G4PreCompoundEmission::rho(const G4double p, const G4double h, const G4double g, 
				    const G4double E, const G4double Ef) const
{
  G4double fact[30];
  fact[0] = 1;
  for (G4int n = 1; n < 21; n++) {
    fact[n] = fact[n-1]*G4double(n);
  }
	
  G4double aph = (p*p + h*h + p - 3.0*h)/(4.0*g);
	
  G4double tot = 0.0;
  for (G4int j = 0; j <= h; j++) {
    G4double t1 = pow(-1.0, G4double(j));
    G4double t2 = fact[j]/ (fact[G4int(h)-j]*fact[G4int(h)]);
    G4double t3 = E - G4double(j)*Ef - aph;
    if (t3 < 0.0) t3 = 0.0;
    t3 = pow(t3,p+h-1);
    tot += t1*t2*t3;
  }
	
  tot *= pow(g,p+h)/(fact[G4int(p)]*fact[G4int(h)]*fact[G4int(p+h)-1]);

  return tot;
}

//  G4double G4PreCompoundEmission::bessi0(const G4double x) const
//    // Returns the modified Bessel function I_0(x) for any real x.
//  {
//    G4double ax,ans; 
//    G4double y;     

//    if ((ax=fabs(x)) < 3.75) {      /* Polynomial fit. */
//      y=x/3.75; 
//      y*=y; 
//      ans=1.0+y*(3.5156229+y*(3.0899424+
//  			    y*(1.2067492+
//  			       y*(0.2659732+
//  				  y*(0.360768e-1+
//  				     y*0.45813e-2))))); 
//    } else {
//      y=3.75/ax; 
//      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+
//  					  y*(0.225319e-2+
//  					     y*(-0.157565e-2+
//  						y*(0.916281e-2+
//  						   y*(-0.2057706e-1+
//  						      y*(0.2635537e-1+
//  							 y*(-0.1647633e-1+
//  							    y*0.392377e-2))))))));
//    } 
//    return ans; 
	
//  }
