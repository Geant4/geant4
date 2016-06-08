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
// $Id: G4CompetitiveFission.cc,v 1.11 2002/01/15 12:27:30 vlara Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)

#include "G4CompetitiveFission.hh"
#include "G4PairingCorrection.hh"

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
}

G4CompetitiveFission::G4CompetitiveFission(const G4CompetitiveFission &right)
{
}

G4CompetitiveFission::~G4CompetitiveFission()
{
    if (MyOwnFissionBarrier) delete theFissionBarrierPtr;

    if (MyOwnFissionProbability) delete theFissionProbabilityPtr;

    if (MyOwnLevelDensity) delete theLevelDensityPtr;
}

const G4CompetitiveFission & G4CompetitiveFission::operator=(const G4CompetitiveFission &right)
{
    G4Exception("G4CompetitiveFission::operator= meant to not be accessable");
    return *this;
}

G4bool G4CompetitiveFission::operator==(const G4CompetitiveFission &right) const
{
    return (this == (G4CompetitiveFission *) &right);
}

G4bool G4CompetitiveFission::operator!=(const G4CompetitiveFission &right) const
{
    return (this != (G4CompetitiveFission *) &right);
}




void G4CompetitiveFission::Initialize(const G4Fragment & fragment)
{
    G4int anA = G4int(fragment.GetA());
    G4int aZ = G4int(fragment.GetZ());
    G4double ExEnergy = fragment.GetExcitationEnergy() - 
	   G4PairingCorrection::GetInstance()->GetFissionPairingCorrection(anA,aZ);
  

    // Saddle point excitation energy ---> A = 65
    // Fission is excluded for A < 65
    if (anA >= 65 && ExEnergy > 0.0) {
	FissionBarrier = theFissionBarrierPtr->FissionBarrier(anA,aZ,ExEnergy);
	MaximalKineticEnergy = ExEnergy - FissionBarrier;
	LevelDensityParameter = theLevelDensityPtr->LevelDensityParameter(anA,aZ,ExEnergy);
	FissionProbability = theFissionProbabilityPtr->EmissionProbability(fragment,MaximalKineticEnergy);
    }
    else {
	MaximalKineticEnergy = -1000.0*MeV;
	LevelDensityParameter = 0.0;
	FissionProbability = 0.0;
    }

    return;
}



G4FragmentVector * G4CompetitiveFission::BreakUp(const G4Fragment & theNucleus)
{
    // Nucleus data
    // Atomic number of nucleus
    G4int A = G4int(theNucleus.GetA());
    // Charge of nucleus
    G4int Z = G4int(theNucleus.GetZ());
    //   Excitation energy (in MeV)
    G4double U = theNucleus.GetExcitationEnergy() - 
	G4PairingCorrection::GetInstance()->GetFissionPairingCorrection(A,Z);
    // Check that U > 0
    if (U <= 0.0) {
	G4FragmentVector * theResult = new  G4FragmentVector;
	theResult->push_back(new G4Fragment(theNucleus));
	return theResult;
    }

    // Atomic Mass of Nucleus (in MeV)
    G4double M = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,A)/MeV;
    // Nucleus Momentum
    G4LorentzVector theNucleusMomentum = theNucleus.GetMomentum();

    // Calculate fission parameters
    G4FissionParameters theParameters(A,Z,U,FissionBarrier);
  
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

    G4int Trials = 0;
    do {

	// First fragment 
	A1 = FissionAtomicNumber(A,theParameters);
	Z1 = FissionCharge(A,Z,A1);
	M1 = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z1,A1);


	// Second Fragment
	A2 = A - A1;
	Z2 = Z - Z1;
	if (A2 < 1 || Z2 < 0) 
	    G4Exception("G4CompetitiveFission::BreakUp: Can't define second fragment! ");
	M2 = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z2,A2)/MeV;

	// Check that fragment masses are less or equal than total energy
	//  if (M1 + M2 > theNucleusMomentum.mag()/MeV)
	if (M1 + M2 > theNucleusMomentum.e()/MeV)
	    G4Exception("G4CompetitiveFission::BreakUp: Fragments Mass > Total Energy");

	// Maximal Kinetic Energy (available energy for fragments)
	//  G4double Tmax = theNucleusMomentum.mag()/MeV - M1 - M2;
	G4double Tmax = M + U - M1 - M2;

	FragmentsKineticEnergy = FissionKineticEnergy( A , Z,
						       A1, Z1,
						       A2, Z2,
						       U , Tmax,
						       theParameters);
    
	// Excitation Energy
	FragmentsExcitationEnergy = Tmax - FragmentsKineticEnergy;
    
    } while (FragmentsExcitationEnergy < 0.0 && Trials++ < 100);
  
  

    if (FragmentsExcitationEnergy <= 0.0) 
	G4Exception("G4CompetitiveFission::BreakItUp: Excitation energy for fragments < 0.0!");
  

  // while (FragmentsExcitationEnergy < 0 && Trials < 100);
  
  // Fragment 1
    G4double U1 = FragmentsExcitationEnergy * (G4double(A1)/G4double(A));
    // Fragment 2
    G4double U2 = FragmentsExcitationEnergy * (G4double(A2)/G4double(A));


    G4double Pmax = sqrt( 2 * ( ( (M1+U1)*(M2+U2) ) /
				( (M1+U1)+(M2+U2) ) ) * FragmentsKineticEnergy);

    G4ParticleMomentum momentum1 = IsotropicVector( Pmax );
    G4ParticleMomentum momentum2( -momentum1 );

    // Perform a Galileo boost for fragments
    momentum1 += (theNucleusMomentum.boostVector() * (M1+U1));
    momentum2 += (theNucleusMomentum.boostVector() * (M2+U2));


    // Create 4-momentum for first fragment
    // Warning!! Energy conservation is broken
    G4LorentzVector FourMomentum1( momentum1 , sqrt(momentum1.mag2() + (M1+U1)*(M1+U1)));

    // Create 4-momentum for second fragment
    // Warning!! Energy conservation is broken
    G4LorentzVector FourMomentum2( momentum2 , sqrt(momentum2.mag2() + (M2+U2)*(M2+U2)));

    // Create Fragments
    G4Fragment * Fragment1 = new G4Fragment( A1, Z1, FourMomentum1);
    if (!Fragment1) G4Exception("G4CompetitiveFission::BreakItUp: Can't create Fragment1! ");
    G4Fragment * Fragment2 = new G4Fragment( A2, Z2, FourMomentum2);
    if (!Fragment2) G4Exception("G4CompetitiveFission::BreakItUp: Can't create Fragment2! ");

#ifdef pctest
    Fragment1->SetCreatorModel(G4String("G4CompetitiveFission"));
    Fragment2->SetCreatorModel(G4String("G4CompetitiveFission"));
#endif
  // Create Fragment Vector
    G4FragmentVector * theResult = new G4FragmentVector;

    theResult->push_back(Fragment1);
    theResult->push_back(Fragment2);

#ifdef debug
    CheckConservation(theNucleus,theResult);
#endif

    return theResult;
}



G4int G4CompetitiveFission::FissionAtomicNumber(const G4int A, const G4FissionParameters & theParam)
    // Calculates the atomic number of a fission product
{

    // For Simplicity reading code
    const G4double A1 = theParam.GetA1();
    const G4double A2 = theParam.GetA2();
    const G4double As = theParam.GetAs();
//    const G4double Sigma1 = theParam.GetSigma1();
    const G4double Sigma2 = theParam.GetSigma2();
    const G4double SigmaS = theParam.GetSigmaS();
    const G4double w = theParam.GetW();

  
//    G4double FasymAsym = 2.0*exp(-((A2-As)*(A2-As))/(2.0*Sigma2*Sigma2)) + 
//	exp(-((A1-As)*(A1-As))/(2.0*Sigma1*Sigma1));

//    G4double FsymA1A2 = exp(-((As-(A1+A2))*(As-(A1+A2)))/(2.0*SigmaS*SigmaS));


    G4double C2A = A2 + 3.72*Sigma2;
    G4double C2S = As + 3.72*SigmaS;
  
    G4double C2 = 0.0;
    if (w > 1000.0 ) C2 = C2S;
    else if (w < 0.001) C2 = C2A;
    else C2 =  G4std::max(C2A,C2S);

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
    G4double m;
    G4double Pm;
    do {
	m = C1+G4UniformRand()*(C2-C1);
	Pm = MassDistribution(m,A,theParam); 
    } while (G4UniformRand() > Pm/MassMax);

    return G4int(m+0.5);
}






G4double G4CompetitiveFission::MassDistribution(const G4double x, const G4double A, 
						const G4FissionParameters & theParam)
    // This method gives mass distribution F(x) = F_{asym}(x)+w*F_{sym}(x)
    // which consist of symmetric and asymmetric sum of gaussians components.
{
    G4double Xsym = exp(-0.5*(x-theParam.GetAs())*(x-theParam.GetAs())/
			(theParam.GetSigmaS()*theParam.GetSigmaS()));

    G4double Xasym = exp(-0.5*(x-theParam.GetA2())*(x-theParam.GetA2())/
			 (theParam.GetSigma2()*theParam.GetSigma2())) + 
	exp(-0.5*(x-(A-theParam.GetA2()))*(x-(A-theParam.GetA2()))/
	    (theParam.GetSigma2()*theParam.GetSigma2())) +
	0.5*exp(-0.5*(x-theParam.GetA1())*(x-theParam.GetA1())/
		(theParam.GetSigma1()*theParam.GetSigma1())) +
	0.5*exp(-0.5*(x-(A-theParam.GetA1()))*(x-(A-theParam.GetA1()))/
		(theParam.GetSigma1()*theParam.GetSigma1()));

    if (theParam.GetW() > 1000) return Xsym;
    else if (theParam.GetW() < 0.001) return Xasym;
    else return theParam.GetW()*Xsym+Xasym;
}




G4int G4CompetitiveFission::FissionCharge(const G4double A,
					  const G4double Z,
					  const G4double Af)
    // Calculates the charge of a fission product for a given atomic number Af
{
    const G4double sigma = 0.6;
    G4double DeltaZ = 0.0;
    if (Af >= 134.0) DeltaZ = -0.45;                    //                      134 <= Af
    else if (A <= (A-134.0)) DeltaZ = 0.45;             // Af <= (A-134) 
    else DeltaZ = -0.45*(Af-(A/2.0))/(134.0-(A/2.0));   //       (A-134) < Af < 134

    G4double Zmean = (Af/A)*Z + DeltaZ;
 
    G4double theZ;
    do {
	theZ = G4RandGauss::shoot(Zmean,sigma);
    } while (theZ  < 1.0 || theZ > (Z-1.0) || theZ > Af);
    //  return static_cast<G4int>(theZ+0.5);
    return G4int(theZ+0.5);
}




G4double G4CompetitiveFission::FissionKineticEnergy(const G4double A, const G4double Z,
						    const G4double Af1, const G4double Zf1,
						    const G4double Af2, const G4double Zf2,
						    const G4double U, const G4double Tmax,
						    const G4FissionParameters & theParam)
    // Gives the kinetic energy of fission products
{
    // Find maximal value of A for fragments
    G4double AfMax = G4std::max(Af1,Af2);
    if (AfMax < (A/2.0)) AfMax = A - AfMax;

    // Weights for symmetric and asymmetric components
    G4double Pas;
    if (theParam.GetW() > 1000) Pas = 0.0;
    else {
	G4double P1 = 0.5*exp(-0.5*(AfMax-theParam.GetA1())*(AfMax-theParam.GetA1())/
			      (theParam.GetSigma1()*theParam.GetSigma1()));

	G4double P2 = exp(-0.5*(AfMax-theParam.GetA2())*(AfMax-theParam.GetA2())/
			  (theParam.GetSigma2()*theParam.GetSigma2()));

	Pas = P1+P2;
    }


    G4double Ps;
    if (theParam.GetW() < 0.001) Ps = 0.0;
    else 
	Ps = theParam.GetW()*exp(-0.5*(AfMax-theParam.GetAs())*(AfMax-theParam.GetAs())/
				 (theParam.GetSigmaS()*theParam.GetSigmaS()));
 


    G4double Psy = Ps/(Pas+Ps);


    // Fission fractions Xsy and Xas formed in symmetric and asymmetric modes
    G4double PPas = theParam.GetSigma1() + 2.0 * theParam.GetSigma2();
    G4double PPsy = theParam.GetW() * theParam.GetSigmaS();
    G4double Xas = PPas / (PPas+PPsy);
    G4double Xsy = PPsy / (PPas+PPsy);


    // Average kinetic energy for symmetric and asymmetric components
    G4double Eaverage = 0.1071*MeV*(Z*Z)/pow(A,1.0/3.0) + 22.2*MeV;


    // Compute maximal average kinetic energy of fragments and Energy Dispersion (sqrt)
    G4double TaverageAfMax;
    G4double ESigma;
    // Select randomly fission mode (symmetric or asymmetric)
    if (G4UniformRand() > Psy) { // Asymmetric Mode
	G4double A11 = theParam.GetA1()-0.7979*theParam.GetSigma1();
	G4double A12 = theParam.GetA1()+0.7979*theParam.GetSigma1();
	G4double A21 = theParam.GetA2()-0.7979*theParam.GetSigma2();
	G4double A22 = theParam.GetA2()+0.7979*theParam.GetSigma2();
	// scale factor
	G4double ScaleFactor = 0.5*theParam.GetSigma1()*(AsymmetricRatio(A,A11)+AsymmetricRatio(A,A12))+
	    theParam.GetSigma2()*(AsymmetricRatio(A,A21)+AsymmetricRatio(A,A22));
	// Compute average kinetic energy for fragment with AfMax
	TaverageAfMax = (Eaverage + 12.5 * Xsy) * (PPas/ScaleFactor) * AsymmetricRatio(A,AfMax);
	ESigma = 10.0*MeV; // MeV

    } else { // Symmetric Mode
	G4double As0 = theParam.GetAs() + 0.7979*theParam.GetSigmaS();
	// scale factor
	G4double ScaleFactor = theParam.GetW()*theParam.GetSigmaS()*SymmetricRatio(A,As0);
	// Compute average kinetic energy for fragment with AfMax
	TaverageAfMax = (Eaverage - 12.5*MeV*Xas) * (PPsy/ScaleFactor) * SymmetricRatio(A,AfMax);
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





G4double G4CompetitiveFission::AsymmetricRatio(const G4double A,const G4double A11)
{
    const G4double B1 = 23.5;
    const G4double A00 = 134.0;
    return Ratio(A,A11,B1,A00);
}

G4double G4CompetitiveFission::SymmetricRatio(const G4double A,const G4double A11)
{
    const G4double B1 = 5.32;
    const G4double A00 = A/2.0;
    return Ratio(A,A11,B1,A00);
}

G4double G4CompetitiveFission::Ratio(const G4double A,const G4double A11,
				     const G4double B1,const G4double A00) 
{
    if (A == 0) G4Exception("G4CompetitiveFission::Ratio: A == 0!");
    if (A11 >= A/2.0 && A11 <= (A00+10.0)) return 1.0-B1*((A11-A00)/A)*((A11-A00)/A);
    else return 1.0-B1*(10.0/A)*(10.0/A)-2.0*(10.0/A)*B1*((A11-A00-10.0)/A);
}





G4ThreeVector G4CompetitiveFission::IsotropicVector(const G4double Magnitude)
    // Samples a isotropic random vectorwith a magnitud given by Magnitude.
    // By default Magnitude = 1.0
{
    G4double CosTheta = 1.0 - 2.0*G4UniformRand();
    G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
    G4double Phi = twopi*G4UniformRand();
    G4ThreeVector Vector(Magnitude*cos(Phi)*SinTheta,
			 Magnitude*sin(Phi)*SinTheta,
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
	ProductsA += G4int((*h)->GetA());
	ProductsZ += G4int((*h)->GetZ());
    }

    if (ProductsA != theInitialState.GetA()) {
	G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4CompetitiveFission.cc: Barionic Number Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial A = " << theInitialState.GetA() 
	       << "   Fragments A = " << ProductsA << "   Diference --> " 
	       << theInitialState.GetA() - ProductsA << G4endl;
    }
    if (ProductsZ != theInitialState.GetZ()) {
	G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4CompetitiveFission.cc: Charge Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial Z = " << theInitialState.GetZ() 
	       << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	       << theInitialState.GetZ() - ProductsZ << G4endl;
    }
    if (abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
	G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4CompetitiveFission.cc: Energy Conservation test for fission fragments" 
	       << G4endl; 
	G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	       << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	       << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
    } 
    if (abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
	abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
	abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
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




