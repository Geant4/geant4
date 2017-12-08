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
// $Id: G4PreCompoundEmission.cc 107062 2017-11-01 15:01:02Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundEmission
//
// Author:         V.Lara
//
// Modified:  
// 15.01.2010 J.M.Quesada  added protection against unphysical values of parameter an 
// 19.01.2010 V.Ivanchenko simplified computation of parameter an, sample cosTheta 
//                         instead of theta; protect all calls to sqrt 
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
//

#include "G4PreCompoundEmission.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4PreCompoundEmissionFactory.hh"
#include "G4HETCEmissionFactory.hh"
#include "G4HadronicException.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

G4PreCompoundEmission::G4PreCompoundEmission()
{
  theFragmentsFactory = new G4PreCompoundEmissionFactory();
  theFragmentsVector = 
    new G4PreCompoundFragmentVector(theFragmentsFactory->GetFragmentVector());
  g4calc = G4Pow::GetInstance();
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters() ;
  fLevelDensity = param->GetLevelDensity();
  fFermiEnergy  = param->GetFermiEnergy();
  fUseAngularGenerator = param->UseAngularGen();
}

G4PreCompoundEmission::~G4PreCompoundEmission()
{
  delete theFragmentsFactory; 
  delete theFragmentsVector; 
}

void G4PreCompoundEmission::SetDefaultModel()
{
  if (theFragmentsFactory) { delete theFragmentsFactory; }
  theFragmentsFactory = new G4PreCompoundEmissionFactory();
  if (theFragmentsVector) {
    theFragmentsVector->SetVector(theFragmentsFactory->GetFragmentVector());
  } else {
    theFragmentsVector = 
      new G4PreCompoundFragmentVector(theFragmentsFactory->GetFragmentVector());
  }
}

void G4PreCompoundEmission::SetHETCModel()
{
  if (theFragmentsFactory) delete theFragmentsFactory;
  theFragmentsFactory = new G4HETCEmissionFactory();
  if (theFragmentsVector) {
    theFragmentsVector->SetVector(theFragmentsFactory->GetFragmentVector());
  } else {
    theFragmentsVector = 
      new G4PreCompoundFragmentVector(theFragmentsFactory->GetFragmentVector());
  }
}

G4ReactionProduct* 
G4PreCompoundEmission::PerformEmission(G4Fragment & aFragment)
{
  // Choose a Fragment for emission
  G4VPreCompoundFragment * thePreFragment = 
    theFragmentsVector->ChooseFragment();
  if (thePreFragment == 0)
    {
      G4cout << "G4PreCompoundEmission::PerformEmission : "
	     << "I couldn't choose a fragment\n"
	     << "while trying to de-excite\n" 
	     << aFragment << G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "");
    }

  //G4cout << "Chosen fragment: " << G4endl;
  //G4cout << *thePreFragment << G4endl;

  // Kinetic Energy of emitted fragment
  G4double kinEnergy = thePreFragment->SampleKineticEnergy(aFragment);
  kinEnergy = std::max(kinEnergy, 0.0);
  
  // Calculate the fragment momentum (three vector)
  if(fUseAngularGenerator) {
    AngularDistribution(thePreFragment,aFragment,kinEnergy);
  } else {
    G4double pmag = 
      std::sqrt(kinEnergy*(kinEnergy + 2.0*thePreFragment->GetNuclearMass()));
    theFinalMomentum = pmag*G4RandomDirection();
  }

  // Mass of emittef fragment
  G4double EmittedMass = thePreFragment->GetNuclearMass();
  // Now we can calculate the four momentum 
  // both options are valid and give the same result but 2nd one is faster
  G4LorentzVector Emitted4Momentum(theFinalMomentum,EmittedMass + kinEnergy);
    
  // Perform Lorentz boost
  G4LorentzVector Rest4Momentum = aFragment.GetMomentum();
  Emitted4Momentum.boost(Rest4Momentum.boostVector());  

  // Set emitted fragment momentum
  thePreFragment->SetMomentum(Emitted4Momentum);	

  // NOW THE RESIDUAL NUCLEUS
  // ------------------------

  Rest4Momentum -= Emitted4Momentum;
    
  // Update nucleus parameters:
  // --------------------------

  // Z and A
  aFragment.SetZandA_asInt(thePreFragment->GetRestZ(),
			   thePreFragment->GetRestA());
    
  // Number of excitons
  aFragment.SetNumberOfParticles(aFragment.GetNumberOfParticles()-
				 thePreFragment->GetA());
  // Number of charges
  aFragment.SetNumberOfCharged(aFragment.GetNumberOfCharged()-
			       thePreFragment->GetZ());
    
  // Update nucleus momentum 
  // A check on consistence of Z, A, and mass will be performed
  aFragment.SetMomentum(Rest4Momentum);
	
  // Create a G4ReactionProduct 
  G4ReactionProduct * MyRP = thePreFragment->GetReactionProduct();

  //  if(kinEnergy < MeV) {
  //  G4cout << "G4PreCompoundEmission::Fragment emitted" << G4endl;
  //  G4cout << *thePreFragment << G4endl;
    // }
  return MyRP;
}

void G4PreCompoundEmission::AngularDistribution(
                            G4VPreCompoundFragment* thePreFragment,
			    const G4Fragment& aFragment,
			    G4double ekin) 
{
  G4int p = aFragment.GetNumberOfParticles();
  G4int h = aFragment.GetNumberOfHoles();
  G4double U = aFragment.GetExcitationEnergy();

  // Emission particle separation energy
  G4double Bemission = thePreFragment->GetBindingEnergy();
	
  G4double gg = (6.0/pi2)*aFragment.GetA_asInt()*fLevelDensity;
	
  // Average exciton energy relative to bottom of nuclear well
  G4double Eav = 2*p*(p+1)/((p+h)*gg);
	
  // Excitation energy relative to the Fermi Level
  G4double Uf = std::max(U - (p - h)*fFermiEnergy , 0.0);
  //  G4double Uf = U - KineticEnergyOfEmittedFragment - Bemission;

  G4double w_num = rho(p+1, h, gg, Uf, fFermiEnergy);
  G4double w_den = rho(p,   h, gg, Uf, fFermiEnergy);
  if (w_num > 0.0 && w_den > 0.0)
    {
      Eav *= (w_num/w_den);
      Eav += - Uf/(p+h) + fFermiEnergy;
    }
  else 
    {
      Eav = fFermiEnergy;
    }
  
  // VI + JMQ 19/01/2010 update computation of the parameter an
  //
  G4double an = 0.0;
  G4double Eeff = ekin + Bemission + fFermiEnergy;
  if(ekin > DBL_MIN && Eeff > DBL_MIN) {

    G4double zeta = std::max(1.0,9.3/std::sqrt(ekin/MeV));
  
    // This should be the projectile energy. If I would know which is 
    // the projectile (proton, neutron) I could remove the binding energy. 
    // But, what happens if INC precedes precompound? This approximation
    // seems to work well enough
    G4double ProjEnergy = aFragment.GetExcitationEnergy();

    an = 3*std::sqrt((ProjEnergy+fFermiEnergy)*Eeff)/(zeta*Eav);

    G4int ne = aFragment.GetNumberOfExcitons() - 1;
    if ( ne > 1 ) { an /= (G4double)ne; }
			
    // protection of exponent
    if ( an > 10. ) { an = 10.; }
  }

  // sample cosine of theta and not theta as in old versions  
  G4double random = G4UniformRand();
  G4double cost;
 
  if(an < 0.1) { cost = 1. - 2*random; }
  else {
    G4double exp2an = G4Exp(-2*an);
    cost = 1. + G4Log(1-random*(1-exp2an))/an;
    if(cost > 1.) { cost = 1.; }
    else if(cost < -1.) {cost = -1.; }
  }  

  G4double phi = CLHEP::twopi*G4UniformRand();
  
  // Calculate the momentum magnitude of emitted fragment 	
  G4double pmag = 
    std::sqrt(ekin*(ekin + 2.0*thePreFragment->GetNuclearMass()));
  
  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  theFinalMomentum.set(pmag*std::cos(phi)*sint,pmag*std::sin(phi)*sint,
		       pmag*cost);

  // theta is the angle wrt the incident direction
  G4ThreeVector theIncidentDirection = aFragment.GetMomentum().vect().unit();
  theFinalMomentum.rotateUz(theIncidentDirection);
}

G4double G4PreCompoundEmission::rho(G4int p, G4int h, G4double gg, 
				    G4double E, G4double Ef) const
{	
  // 25.02.2010 V.Ivanchenko added more protections
  G4double Aph   = (p*p + h*h + p - 3.0*h)/(4.0*gg);
  //  G4double alpha = (p*p + h*h)/(2.0*gg);
  
  if ( E - Aph < 0.0) { return 0.0; }
  
  G4double logConst =  (p+h)*G4Log(gg) 
    - g4calc->logfactorial(p+h-1) - g4calc->logfactorial(p) 
    - g4calc->logfactorial(h);

  // initialise values using j=0

  G4double t1=1;
  G4double t2=1;
  G4double logt3 = (p+h-1) * G4Log(E-Aph) + logConst;
  const G4double logmax = 200.;
  if(logt3 > logmax) { logt3 = logmax; }
  G4double tot = G4Exp( logt3 );

  // and now sum rest of terms
  // 25.02.2010 V.Ivanchenko change while to for loop and cleanup 
  G4double Eeff = E - Aph; 
  for(G4int j=1; j<=h; ++j) 
    {
      Eeff -= Ef;
      if(Eeff < 0.0) { break; }
      t1 *= -1.;
      t2 *= (G4double)(h+1-j)/(G4double)j;
      logt3 = (p+h-1) * G4Log( Eeff) + logConst;
      if(logt3 > logmax) { logt3 = logmax; }
      tot += t1*t2*G4Exp(logt3);
    }
        
  return tot;
}
