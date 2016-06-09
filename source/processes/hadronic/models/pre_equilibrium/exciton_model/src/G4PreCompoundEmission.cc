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
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 


#include "G4PreCompoundEmission.hh"
#include "G4PreCompoundParameters.hh"

#include "G4PreCompoundEmissionFactory.hh"
#include "G4HETCEmissionFactory.hh"

const G4PreCompoundEmission & G4PreCompoundEmission::operator=(const G4PreCompoundEmission &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundEmission::operator= meant to not be accessable");
  return *this;
}


  G4bool G4PreCompoundEmission::operator==(const G4PreCompoundEmission &) const
{
  return false;
}

    G4bool G4PreCompoundEmission::operator!=(const G4PreCompoundEmission &) const
{
  return true;
}

G4PreCompoundEmission::G4PreCompoundEmission()
{
  theFragmentsFactory = new G4PreCompoundEmissionFactory();
  theFragmentsVector = new G4PreCompoundFragmentVector(theFragmentsFactory->GetFragmentVector());
}

G4PreCompoundEmission::~G4PreCompoundEmission()
{
  if (theFragmentsFactory) delete theFragmentsFactory;
  if (theFragmentsVector) delete theFragmentsVector;
}

void G4PreCompoundEmission::SetDefaultModel()
{
  if (theFragmentsFactory) delete theFragmentsFactory;
  theFragmentsFactory = new G4PreCompoundEmissionFactory();
  if (theFragmentsVector) 
    {
      theFragmentsVector->SetVector(theFragmentsFactory->GetFragmentVector());
    }
  else 
    {
      theFragmentsVector = new G4PreCompoundFragmentVector(theFragmentsFactory->GetFragmentVector());
    }
  theFragmentsVector->ResetStage();
  return;
}

void G4PreCompoundEmission::SetHETCModel()
{
  if (theFragmentsFactory) delete theFragmentsFactory;
  theFragmentsFactory = new G4HETCEmissionFactory();
  if (theFragmentsVector) 
    {
      theFragmentsVector->SetVector(theFragmentsFactory->GetFragmentVector());
    }
  else 
    {
      theFragmentsVector = new G4PreCompoundFragmentVector(theFragmentsFactory->GetFragmentVector());
    }
  theFragmentsVector->ResetStage();
  return;
}



G4ReactionProduct * G4PreCompoundEmission::PerformEmission(G4Fragment & aFragment)
{
#ifdef debug
  G4Fragment InitialState(aFragment);
#endif
  // Choose a Fragment for emission
  G4VPreCompoundFragment * theFragment = theFragmentsVector->ChooseFragment();
  if (theFragment == 0)
    {
      G4cerr <<  "G4PreCompoundEmission::PerformEmission : I couldn't choose a fragment\n"
	     << "while trying to de-excite\n" 
	     << aFragment << '\n';
      throw G4HadronicException(__FILE__, __LINE__, "");
    }
  // Kinetic Energy of emitted fragment
  G4double KineticEnergyOfEmittedFragment = theFragment->GetKineticEnergy(aFragment);
  
  // Calculate the fragment momentum (three vector)
  G4ThreeVector momentum = AngularDistribution(theFragment,aFragment,KineticEnergyOfEmittedFragment);
  
  // Mass of emittef fragment
  G4double EmittedMass = theFragment->GetNuclearMass();
  
  // Now we can calculate the four momentum 
  // both options are valid and give the same result but 2nd one is faster
  // G4LorentzVector EmittedMomentum(momentum,std::sqrt(momentum.mag2()+EmittedMass*EmittedMass));
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
    //    aFragment.GetGroundStateMass() + aFragment.GetExcitationEnergy() // initial energy in cm
    aFragment.GetMomentum().m()
    - (EmittedMass+KineticEnergyOfEmittedFragment); 

  // Then the four momentum for residual is 
  G4LorentzVector RestMomentum(-momentum,ResidualEcm);
  // This could save a Lorentz boost
  // G4LorentzVector RestMomentum2(aFragment.GetMomentum()-EmittedMomentum);

  // Just for test
  // Excitation energy
  //  G4double anU = ResidualEcm - theFragment->GetRestNuclearMass();
  // This is equivalent
  //  G4double anU = theFragment->GetMaximalKineticEnergy() - KineticEnergyOfEmittedFragment + 
  //    theFragment->GetCoulombBarrier();
    
  // check that Excitation energy is >= 0
  G4double anU = RestMomentum.m()-theFragment->GetRestNuclearMass();
  if (anU < 0.0) throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundModel::DeExcite: Excitation energy less than 0!");
    
    
    
  // Update nucleus parameters:
  // --------------------------

  // Number of excitons
  aFragment.SetNumberOfParticles(aFragment.GetNumberOfParticles()-
				 static_cast<G4int>(theFragment->GetA()));
  // Number of charges
  aFragment.SetNumberOfCharged(aFragment.GetNumberOfCharged()-
			       static_cast<G4int>(theFragment->GetZ()));
    
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
#ifdef PRECOMPOUND_TEST
  MyRP->SetCreatorModel("G4PreCompoundModel");
#endif
#ifdef debug
  CheckConservation(InitialState,aFragment,MyRP);
#endif
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
  //  G4EvaporationLevelDensityParameter theLDP;
  //  G4double g = (6.0/pi2)*aFragment.GetA()*
  //    theLDP.LevelDensityParameter(static_cast<G4int>(aFragment.GetA()),static_cast<G4int>(aFragment.GetZ()),U);
  G4double g = (6.0/pi2)*aFragment.GetA()*
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
	
  // Average exciton energy relative to bottom of nuclear well
  G4double Eav = 2.0*p*(p+1.0)/((p+h)*g);
	
  // Excitation energy relative to the Fermi Level
  G4double Uf = std::max(U - (p - h)*Ef , 0.0);
  //  G4double Uf = U - KineticEnergyOfEmittedFragment - Bemission;



  G4double w_num = rho(p+1,h,g,Uf,Ef);
  G4double w_den = rho(p,h,g,Uf,Ef);
  if (w_num > 0.0 && w_den > 0.0)
    {
      Eav *= (w_num/w_den);
      Eav += - Uf/(p+h) + Ef;
    }
  else 
    {
      Eav = Ef;
    }
  
  G4double zeta = std::max(1.0,9.3/std::sqrt(KineticEnergyOfEmittedFragment/MeV));
	
  G4double an = 3.0*std::sqrt((ProjEnergy+Ef)*(KineticEnergyOfEmittedFragment+Bemission+Ef));
  if (aFragment.GetNumberOfExcitons() == 1)
    {
      an /= (zeta*2.0*aFragment.GetNumberOfExcitons()*Eav);
    }
  else
    {
      an /= (zeta*(aFragment.GetNumberOfExcitons()-1.0)*Eav);
    }
			
			
  G4double expan = std::exp(an);
	
  G4double theta = std::acos(std::log(expan-G4UniformRand()*(expan-1.0/expan))/an);
	
  G4double phi = twopi*G4UniformRand();
  
  // Calculate the momentum magnitude of emitted fragment 	
  G4double EmittedMass = theFragment->GetNuclearMass();
  G4double pmag = std::sqrt(KineticEnergyOfEmittedFragment*(KineticEnergyOfEmittedFragment+2.0*EmittedMass));
  
  
  G4double sinTheta = std::sin(theta);
  //  G4double cosTheta = std::sqrt(1.0-sinTheta*sinTheta);
  G4double cosTheta = std::cos(theta);

  G4ThreeVector momentum(pmag*std::cos(phi)*sinTheta,pmag*std::sin(phi)*sinTheta,pmag*cosTheta);
  // theta is the angle wrt the incident direction
  momentum.rotateUz(theIncidentDirection);

  return momentum;
}

G4double G4PreCompoundEmission::rho(const G4double p, const G4double h, const G4double g, 
				    const G4double E, const G4double Ef) const
{
	
  G4double Aph = (p*p + h*h + p - 3.0*h)/(4.0*g);
  G4double alpha = (p*p+h*h)/(2.0*g);
  
  if ( (E-alpha) < 0 ) return 0;

  G4double factp=factorial(p);

  G4double facth=factorial(h);

  G4double factph=factorial(p+h-1);
  
  G4double logConst =  (p+h)*std::log(g) - std::log (factph) - std::log(factp) - std::log(facth);

// initialise values using j=0

  G4double t1=1;
  G4double t2=1;
  G4double logt3=(p+h-1) * std::log(E-Aph);
  G4double tot = std::exp( logt3 + logConst );

// and now sum rest of terms 
  G4int j(1);  
  while ( (j <= h) && ((E - alpha - j*Ef) > 0.0) ) 
    {
	  t1 *= -1.;
	  t2 *= (h+1-j)/j;
	  logt3 = (p+h-1) * std::log( E - j*Ef - Aph) + logConst;
	  G4double t3 = std::exp(logt3);
	  tot += t1*t2*t3;
	  j++;
    }
        
  return tot;
}



G4double G4PreCompoundEmission::factorial(G4double a) const
{
  // Values of factorial function from 0 to 60
  const G4int factablesize = 61;
  static const G4double fact[factablesize] = 
    {
      1.0, // 0!
      1.0, // 1!
      2.0, // 2!
      6.0, // 3!
      24.0, // 4!
      120.0, // 5!
      720.0, // 6!
      5040.0, // 7!
      40320.0, // 8!
      362880.0, // 9!
      3628800.0, // 10!
      39916800.0, // 11!
      479001600.0, // 12!
      6227020800.0, // 13!
      87178291200.0, // 14!
      1307674368000.0, // 15!
      20922789888000.0, // 16!
      355687428096000.0, // 17!
      6402373705728000.0, // 18!
      121645100408832000.0, // 19!
      2432902008176640000.0, // 20!
      51090942171709440000.0, // 21!
      1124000727777607680000.0, // 22!
      25852016738884976640000.0, // 23!
      620448401733239439360000.0, // 24!
      15511210043330985984000000.0, // 25!
      403291461126605635584000000.0, // 26!
      10888869450418352160768000000.0, // 27!
      304888344611713860501504000000.0, // 28!
      8841761993739701954543616000000.0, // 29!
      265252859812191058636308480000000.0, // 30!
      8222838654177922817725562880000000.0, // 31!
      263130836933693530167218012160000000.0, // 32!
      8683317618811886495518194401280000000.0, // 33!
      295232799039604140847618609643520000000.0, // 34!
      10333147966386144929666651337523200000000.0, // 35!
      371993326789901217467999448150835200000000.0, // 36!
      13763753091226345046315979581580902400000000.0, // 37!
      523022617466601111760007224100074291200000000.0, // 38!
      20397882081197443358640281739902897356800000000.0, // 39!
      815915283247897734345611269596115894272000000000.0, // 40!
      33452526613163807108170062053440751665152000000000.0, // 41!
      1405006117752879898543142606244511569936384000000000.0, // 42!
      60415263063373835637355132068513997507264512000000000.0, // 43!
      2658271574788448768043625811014615890319638528000000000.0, // 44!
      119622220865480194561963161495657715064383733760000000000.0, // 45!
      5502622159812088949850305428800254892961651752960000000000.0, // 46!
      258623241511168180642964355153611979969197632389120000000000.0, // 47!
      12413915592536072670862289047373375038521486354677760000000000.0, // 48!
      608281864034267560872252163321295376887552831379210240000000000.0, // 49!
      30414093201713378043612608166064768844377641568960512000000000000.0, // 50!
      1551118753287382280224243016469303211063259720016986112000000000000.0, // 51!
      80658175170943878571660636856403766975289505440883277824000000000000.0, // 52!
      4274883284060025564298013753389399649690343788366813724672000000000000.0, // 53!
      230843697339241380472092742683027581083278564571807941132288000000000000.0, // 54!
      12696403353658275925965100847566516959580321051449436762275840000000000000.0, // 55!
      710998587804863451854045647463724949736497978881168458687447040000000000000.0, // 56!
      40526919504877216755680601905432322134980384796226602145184481280000000000000.0, // 57!
      2350561331282878571829474910515074683828862318181142924420699914240000000000000.0, // 58!
      138683118545689835737939019720389406345902876772687432540821294940160000000000000.0, // 59!
      8320987112741390144276341183223364380754172606361245952449277696409600000000000000.0  // 60!
    };
  //    fact[0] = 1;
  //    for (G4int n = 1; n < 21; n++) {
  //      fact[n] = fact[n-1]*static_cast<G4double>(n); 
  //    }
  G4double result(0.0);
  G4int ia = static_cast<G4int>(a);
  if (ia < factablesize) 
    {
      result = fact[ia];
    }
  else
    {
      result = fact[factablesize-1];
      for (G4int n = factablesize; n < ia+1; ++n)
        {
          result *= static_cast<G4double>(n);
        }
    }
    
    return result;
}


#ifdef debug
void G4PreCompoundEmission::CheckConservation(const G4Fragment & theInitialState,
					      const G4Fragment & theResidual,
					      G4ReactionProduct * theEmitted) const
{
  G4double ProductsEnergy = theEmitted->GetTotalEnergy() + theResidual.GetMomentum().e();
  G4ThreeVector ProductsMomentum(theEmitted->GetMomentum()+theResidual.GetMomentum().vect());
  G4int ProductsA = theEmitted->GetDefinition()->GetBaryonNumber() + theResidual.GetA();
  G4int ProductsZ = theEmitted->GetDefinition()->GetPDGCharge() + theResidual.GetZ();

  if (ProductsA != theInitialState.GetA()) {
    G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundEmission.cc: Barionic Number Conservation"
	   << G4endl; 
    G4cout << "Initial A = " << theInitialState.GetA() 
	   << "   Fragments A = " << ProductsA << "   Diference --> " 
	   << theInitialState.GetA() - ProductsA << G4endl;
  }
  if (ProductsZ != theInitialState.GetZ()) {
    G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundEmission.cc: Charge Conservation test"
	   << G4endl; 
    G4cout << "Initial Z = " << theInitialState.GetZ() 
	   << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	   << theInitialState.GetZ() - ProductsZ << G4endl;
  }
  if (std::abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 10.0*eV) {
    G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundEmission.cc: Energy Conservation test" 
	   << G4endl; 
    G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	   << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	   << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
  } 
  if (std::abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 10.0*eV || 
      std::abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 10.0*eV ||
      std::abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 10.0*eV) {
    G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundEmission.cc: Momentum Conservation test" 
	   << G4endl; 
    G4cout << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	   << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	   << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV" << G4endl;
  }
  return;
}

#endif
