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
// $Id: G4CompetitiveFission.cc 85841 2014-11-05 15:35:06Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// J. M. Quesada (March 2009). Bugs fixed:
//          - Full relativistic calculation (Lorentz boosts)
//          - Fission pairing energy is included in fragment excitation energies
// Now Energy and momentum are conserved in fission 

#include "G4CompetitiveFission.hh"
#include "G4PairingCorrection.hh"
#include "G4ParticleMomentum.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4CompetitiveFission::G4CompetitiveFission() : G4VEvaporationChannel("fission")
{
  theFissionBarrierPtr = new G4FissionBarrier;
  MyOwnFissionBarrier = true;

  theFissionProbabilityPtr = new G4FissionProbability;
  MyOwnFissionProbability = true;
  
  theLevelDensityPtr = new G4FissionLevelDensityParameter;
  MyOwnLevelDensity = true;

  MaximalKineticEnergy = -1000.0*MeV;
  FissionBarrier = 0.0;
  FissionProbability = 0.0;
  LevelDensityParameter = 0.0;
  pairingCorrection = G4PairingCorrection::GetInstance();
}

G4CompetitiveFission::~G4CompetitiveFission()
{
  if (MyOwnFissionBarrier) delete theFissionBarrierPtr;
  if (MyOwnFissionProbability) delete theFissionProbabilityPtr;
  if (MyOwnLevelDensity) delete theLevelDensityPtr;
}

G4double G4CompetitiveFission::GetEmissionProbability(G4Fragment* fragment)
{
  G4int anA = fragment->GetA_asInt();
  G4int aZ  = fragment->GetZ_asInt();
  G4double ExEnergy = fragment->GetExcitationEnergy() - 
    pairingCorrection->GetFissionPairingCorrection(anA,aZ);
  

  // Saddle point excitation energy ---> A = 65
  // Fission is excluded for A < 65
  if (anA >= 65 && ExEnergy > 0.0) {
    FissionBarrier = theFissionBarrierPtr->FissionBarrier(anA,aZ,ExEnergy);
    MaximalKineticEnergy = ExEnergy - FissionBarrier;
    LevelDensityParameter = 
      theLevelDensityPtr->LevelDensityParameter(anA,aZ,ExEnergy);
    FissionProbability = 
      theFissionProbabilityPtr->EmissionProbability(*fragment,MaximalKineticEnergy);
    }
  else {
    MaximalKineticEnergy = -1000.0*MeV;
    LevelDensityParameter = 0.0;
    FissionProbability = 0.0;
  }
  return FissionProbability;
}

G4FragmentVector * G4CompetitiveFission::BreakUp(const G4Fragment & theNucleus)
{
  G4FragmentVector * theResult = new G4FragmentVector();
  G4Fragment* frag0 = new G4Fragment(theNucleus);
  G4Fragment* frag1 = EmittedFragment(frag0);
  if(frag1) { theResult->push_back(frag1); }
  theResult->push_back(frag0);
  return theResult;
}

G4Fragment* G4CompetitiveFission::EmittedFragment(G4Fragment* theNucleus)
{
  G4Fragment * Fragment1 = 0; 
  // Nucleus data
  // Atomic number of nucleus
  G4int A = theNucleus->GetA_asInt();
  // Charge of nucleus
  G4int Z = theNucleus->GetZ_asInt();
  //   Excitation energy (in MeV)
  G4double U = theNucleus->GetExcitationEnergy();
  G4double pcorr = pairingCorrection->GetFissionPairingCorrection(A,Z);
  if (U <= pcorr) { return Fragment1; }

  // Atomic Mass of Nucleus (in MeV)
  G4double M = theNucleus->GetGroundStateMass();

  // Nucleus Momentum
  G4LorentzVector theNucleusMomentum = theNucleus->GetMomentum();

  // Calculate fission parameters
  G4FissionParameters theParameters(A,Z,U-pcorr,FissionBarrier);
  
  // First fragment
  G4int A1 = 0;
  G4int Z1 = 0;
  G4double M1 = 0.0;

  // Second fragment
  G4int A2 = 0;
  G4int Z2 = 0;
  G4double M2 = 0.0;

  G4double FragmentsExcitationEnergy = 0.0;
  G4double FragmentsKineticEnergy = 0.0;

  //JMQ 04/03/09 It will be used latter to fix the bug in energy conservation
  G4double FissionPairingEnergy=
    pairingCorrection->GetFissionPairingCorrection(A,Z);

  G4int Trials = 0;
  do {

    // First fragment 
    A1 = FissionAtomicNumber(A,theParameters);
    Z1 = FissionCharge(A,Z,A1);
    M1 = G4NucleiProperties::GetNuclearMass(A1, Z1);

    // Second Fragment
    A2 = A - A1;
    Z2 = Z - Z1;
    if (A2 < 1 || Z2 < 0 || Z2 > A2) {
      FragmentsExcitationEnergy = -1.0;
      continue;
    }
    M2 = G4NucleiProperties::GetNuclearMass(A2, Z2);
    // Maximal Kinetic Energy (available energy for fragments)
    G4double Tmax = M + U - M1 - M2;

    // Check that fragment masses are less or equal than total energy
    if (Tmax < 0.0) {
      FragmentsExcitationEnergy = -1.0;
      continue;
    }

    FragmentsKineticEnergy = FissionKineticEnergy( A , Z,
						   A1, Z1,
						   A2, Z2,
						   U , Tmax,
						   theParameters);
    
    // Excitation Energy
    //	FragmentsExcitationEnergy = Tmax - FragmentsKineticEnergy;
    // JMQ 04/03/09 BUG FIXED: in order to fulfill energy conservation the
    // fragments carry the fission pairing energy in form of 
    // excitation energy

    FragmentsExcitationEnergy = 
      Tmax - FragmentsKineticEnergy+FissionPairingEnergy;

  } while (FragmentsExcitationEnergy < 0.0 && Trials++ < 100);
    
  if (FragmentsExcitationEnergy <= 0.0) { 
    throw G4HadronicException(__FILE__, __LINE__, 
      "G4CompetitiveFission::BreakItUp: Excitation energy for fragments < 0.0!");
  }

  // Fragment 1
  M1 += FragmentsExcitationEnergy * A1/static_cast<G4double>(A);
  // Fragment 2
  M2 += FragmentsExcitationEnergy * A2/static_cast<G4double>(A);
  // primary
  M += U;

  G4double etot1 = (M*M - M2*M2 + M1*M1)/(2*M);
  G4ParticleMomentum Momentum1(IsotropicVector(std::sqrt((etot1 - M1)*(etot1+M1))));
  G4LorentzVector FourMomentum1(Momentum1, etot1);
  FourMomentum1.boost(theNucleusMomentum.boostVector());
    
  // Create Fragments
  Fragment1 = new G4Fragment( A1, Z1, FourMomentum1);
  theNucleusMomentum -= FourMomentum1;
  theNucleus->SetZandA_asInt(Z2, A2);
  theNucleus->SetMomentum(theNucleusMomentum);
  return Fragment1;
}

G4int 
G4CompetitiveFission::FissionAtomicNumber(G4int A, 
					  const G4FissionParameters & theParam)
  // Calculates the atomic number of a fission product
{

  // For Simplicity reading code
  G4double A1 = theParam.GetA1();
  G4double A2 = theParam.GetA2();
  G4double As = theParam.GetAs();
  G4double Sigma2 = theParam.GetSigma2();
  G4double SigmaS = theParam.GetSigmaS();
  G4double w = theParam.GetW();
  
  //    G4double FasymAsym = 2.0*std::exp(-((A2-As)*(A2-As))/(2.0*Sigma2*Sigma2)) + 
  //	std::exp(-((A1-As)*(A1-As))/(2.0*Sigma1*Sigma1));

  //    G4double FsymA1A2 = std::exp(-((As-(A1+A2))*(As-(A1+A2)))/(2.0*SigmaS*SigmaS));

  G4double C2A = A2 + 3.72*Sigma2;
  G4double C2S = As + 3.72*SigmaS;
  
  G4double C2 = 0.0;
  if (w > 1000.0 ) C2 = C2S;
  else if (w < 0.001) C2 = C2A;
  else C2 =  std::max(C2A,C2S);

  G4double C1 = A-C2;
  if (C1 < 30.0) {
    C2 = A-30.0;
    C1 = 30.0;
  }

  G4double Am1 = (As + A1)/2.0;
  G4double Am2 = (A1 + A2)/2.0;

  // Get Mass distributions as sum of symmetric and asymmetric Gasussians
  G4double Mass1 = MassDistribution(As,A,theParam); 
  G4double Mass2 = MassDistribution(Am1,A,theParam); 
  G4double Mass3 = MassDistribution(A1,A,theParam); 
  G4double Mass4 = MassDistribution(Am2,A,theParam); 
  G4double Mass5 = MassDistribution(A2,A,theParam); 
  // get maximal value among Mass1,...,Mass5
  G4double MassMax = Mass1;
  if (Mass2 > MassMax) MassMax = Mass2;
  if (Mass3 > MassMax) MassMax = Mass3;
  if (Mass4 > MassMax) MassMax = Mass4;
  if (Mass5 > MassMax) MassMax = Mass5;

  // Sample a fragment mass number, which lies between C1 and C2
  G4double xm;
  G4double Pm;
  do {
    xm = C1+G4UniformRand()*(C2-C1);
    Pm = MassDistribution(xm,A,theParam); 
  } while (MassMax*G4UniformRand() > Pm);
  G4int ires = G4lrint(xm);

  return ires;
}

G4double 
G4CompetitiveFission::MassDistribution(G4double x, G4double A, 
				       const G4FissionParameters & theParam)
  // This method gives mass distribution F(x) = F_{asym}(x)+w*F_{sym}(x)
  // which consist of symmetric and asymmetric sum of gaussians components.
{
  G4double Xsym = std::exp(-0.5*(x-theParam.GetAs())*(x-theParam.GetAs())/
			   (theParam.GetSigmaS()*theParam.GetSigmaS()));

  G4double Xasym = std::exp(-0.5*(x-theParam.GetA2())*(x-theParam.GetA2())/
			    (theParam.GetSigma2()*theParam.GetSigma2())) + 
    std::exp(-0.5*(x-(A-theParam.GetA2()))*(x-(A-theParam.GetA2()))/
	     (theParam.GetSigma2()*theParam.GetSigma2())) +
    0.5*std::exp(-0.5*(x-theParam.GetA1())*(x-theParam.GetA1())/
		 (theParam.GetSigma1()*theParam.GetSigma1())) +
    0.5*std::exp(-0.5*(x-(A-theParam.GetA1()))*(x-(A-theParam.GetA1()))/
		 (theParam.GetSigma1()*theParam.GetSigma1()));

  if (theParam.GetW() > 1000) return Xsym;
  else if (theParam.GetW() < 0.001) return Xasym;
  else return theParam.GetW()*Xsym+Xasym;
}

G4int G4CompetitiveFission::FissionCharge(G4double A, G4double Z,
					  G4double Af)
  // Calculates the charge of a fission product for a given atomic number Af
{
  static const G4double sigma = 0.6;
  G4double DeltaZ = 0.0;
  if (Af >= 134.0) DeltaZ = -0.45;                    //                      134 <= Af
  else if (Af <= (A-134.0)) DeltaZ = 0.45;             // Af <= (A-134) 
  else DeltaZ = -0.45*(Af-(A/2.0))/(134.0-(A/2.0));   //       (A-134) < Af < 134

  G4double Zmean = (Af/A)*Z + DeltaZ;
 
  G4double theZ;
  do {
    theZ = G4RandGauss::shoot(Zmean,sigma);
  } while (theZ  < 1.0 || theZ > (Z-1.0) || theZ > Af);
  //  return static_cast<G4int>(theZ+0.5);
  return static_cast<G4int>(theZ+0.5);
}

G4double 
G4CompetitiveFission::FissionKineticEnergy(G4int A, G4int Z,
					   G4double Af1, G4double /*Zf1*/,
					   G4double Af2, G4double /*Zf2*/,
					   G4double /*U*/, G4double Tmax,
					   const G4FissionParameters & theParam)
  // Gives the kinetic energy of fission products
{
  // Find maximal value of A for fragments
  G4double AfMax = std::max(Af1,Af2);
  if (AfMax < (A/2.0)) AfMax = A - AfMax;

  // Weights for symmetric and asymmetric components
  G4double Pas;
  if (theParam.GetW() > 1000) Pas = 0.0;
  else {
    G4double P1 = 0.5*std::exp(-0.5*(AfMax-theParam.GetA1())*(AfMax-theParam.GetA1())/
			       (theParam.GetSigma1()*theParam.GetSigma1()));

    G4double P2 = std::exp(-0.5*(AfMax-theParam.GetA2())*(AfMax-theParam.GetA2())/
			   (theParam.GetSigma2()*theParam.GetSigma2()));

    Pas = P1+P2;
  }

  G4double Ps;
  if (theParam.GetW() < 0.001) Ps = 0.0;
  else {
    Ps = theParam.GetW()*std::exp(-0.5*(AfMax-theParam.GetAs())*(AfMax-theParam.GetAs())/
				  (theParam.GetSigmaS()*theParam.GetSigmaS()));
  }
  G4double Psy = Ps/(Pas+Ps);

  // Fission fractions Xsy and Xas formed in symmetric and asymmetric modes
  G4double PPas = theParam.GetSigma1() + 2.0 * theParam.GetSigma2();
  G4double PPsy = theParam.GetW() * theParam.GetSigmaS();
  G4double Xas = PPas / (PPas+PPsy);
  G4double Xsy = PPsy / (PPas+PPsy);

  // Average kinetic energy for symmetric and asymmetric components
  G4double Eaverage = 0.1071*MeV*(Z*Z)/G4Pow::GetInstance()->Z13(A) + 22.2*MeV;

  // Compute maximal average kinetic energy of fragments and Energy Dispersion (sqrt)
  G4double TaverageAfMax;
  G4double ESigma = 10*MeV;
  // Select randomly fission mode (symmetric or asymmetric)
  if (G4UniformRand() > Psy) { // Asymmetric Mode
    G4double A11 = theParam.GetA1()-0.7979*theParam.GetSigma1();
    G4double A12 = theParam.GetA1()+0.7979*theParam.GetSigma1();
    G4double A21 = theParam.GetA2()-0.7979*theParam.GetSigma2();
    G4double A22 = theParam.GetA2()+0.7979*theParam.GetSigma2();
    // scale factor
    G4double ScaleFactor = 0.5*theParam.GetSigma1()*
      (AsymmetricRatio(A,A11)+AsymmetricRatio(A,A12))+
      theParam.GetSigma2()*(AsymmetricRatio(A,A21)+AsymmetricRatio(A,A22));
    // Compute average kinetic energy for fragment with AfMax
    TaverageAfMax = (Eaverage + 12.5 * Xsy) * (PPas/ScaleFactor) * 
      AsymmetricRatio(A,AfMax);

  } else { // Symmetric Mode
    G4double As0 = theParam.GetAs() + 0.7979*theParam.GetSigmaS();
    // scale factor
    G4double ScaleFactor = theParam.GetW()*theParam.GetSigmaS()*SymmetricRatio(A,As0);
    // Compute average kinetic energy for fragment with AfMax
    TaverageAfMax = (Eaverage - 12.5*MeV*Xas) * (PPsy/ScaleFactor) * 
      SymmetricRatio(A,AfMax);
    ESigma = 8.0*MeV;
  }


  // Select randomly, in accordance with Gaussian distribution, fragment kinetic energy
  G4double KineticEnergy;
  G4int i = 0;
  do {
    KineticEnergy = G4RandGauss::shoot(TaverageAfMax,ESigma);
    if (i++ > 100) return Eaverage;
  } while (KineticEnergy < Eaverage-3.72*ESigma || 
	   KineticEnergy > Eaverage+3.72*ESigma ||
	   KineticEnergy > Tmax);
  
  return KineticEnergy;
}

G4double G4CompetitiveFission::AsymmetricRatio(G4int A, G4double A11)
{
  static const G4double B1 = 23.5;
  static const G4double A00 = 134.0;
  return Ratio(G4double(A),A11,B1,A00);
}

G4double G4CompetitiveFission::SymmetricRatio(G4int A, G4double A11)
{
  static const G4double B1 = 5.32;
  const G4double A00 = A/2.0;
  return Ratio(G4double(A),A11,B1,A00);
}

G4double G4CompetitiveFission::Ratio(G4double A, G4double A11,
				     G4double B1, G4double A00) 
{
  if (A == 0.0) {
    throw G4HadronicException(__FILE__, __LINE__, 
			      "G4CompetitiveFission::Ratio: A == 0!");
  }
  if (A11 >= A/2.0 && A11 <= (A00+10.0)) {
    return 1.0-B1*((A11-A00)/A)*((A11-A00)/A);
  } else {
    return 1.0-B1*(10.0/A)*(10.0/A)-2.0*(10.0/A)*B1*((A11-A00-10.0)/A);
  }
}

G4ThreeVector G4CompetitiveFission::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = std::sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*std::cos(Phi)*SinTheta,
		       Magnitude*std::sin(Phi)*SinTheta,
		       Magnitude*CosTheta);
  return Vector;
}

#ifdef debug
void G4CompetitiveFission::CheckConservation(const G4Fragment & theInitialState,
					     G4FragmentVector * Result) const
{
    G4double ProductsEnergy =0;
    G4ThreeVector ProductsMomentum;
    G4int ProductsA = 0;
    G4int ProductsZ = 0;
    G4FragmentVector::iterator h;
    for (h = Result->begin(); h != Result->end(); h++) {
	G4LorentzVector tmp = (*h)->GetMomentum();
	ProductsEnergy += tmp.e();
	ProductsMomentum += tmp.vect();
	ProductsA += (*h)->GetA_asInt();
	ProductsZ += (*h)->GetZ_asInt();
    }

    if (ProductsA != theInitialState.GetA_asInt()) {
	G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" 
	       << G4endl;
	G4cout << "G4CompetitiveFission: Baryon Number Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial A = " << theInitialState.GetA_asInt() 
	       << "   Fragments A = " << ProductsA << "   Diference --> " 
	       << theInitialState.GetA_asInt() - ProductsA << G4endl;
    }
    if (ProductsZ != theInitialState.GetZ_asInt()) {
	G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4CompetitiveFission.cc: Charge Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial Z = " << theInitialState.GetZ_asInt() 
	       << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	       << theInitialState.GetZ() - ProductsZ << G4endl;
    }
    if (std::fabs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
	G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4CompetitiveFission.cc: Energy Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	       << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	       << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
    } 
    if (std::fabs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
	std::fabs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
	std::fabs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
	G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4CompetitiveFission.cc: Momentum Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	       << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	       << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV" << G4endl;
    }
    return;
}
#endif




