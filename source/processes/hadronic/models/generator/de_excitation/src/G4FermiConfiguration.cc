// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

// V. Lara (Apr 1999)
// Corrected a bug in calculation of probabilities found by N. Ameline
// 


#include "G4FermiConfiguration.hh"

// Kappa = V/V_0 it is used in calculation of Coulomb energy
// Kappa is adimensional
const G4double G4FermiConfiguration::Kappa = 1.0;

// r0 is the nuclear radius
const G4double G4FermiConfiguration::r0 = 1.3*fermi;

//                                                       A  Z  Pol  ExcitE
G4StableFermiFragment G4FermiConfiguration::Fragment00(  1, 0,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment01(  1, 1,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment02(  2, 1,  3,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment03(  3, 1,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment04(  3, 2,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment05(  4, 2,  1,  0.00*keV );
G4He5FermiFragment    G4FermiConfiguration::Fragment06(  5, 2,  4, 16.76*keV ); // He5
G4Li5FermiFragment    G4FermiConfiguration::Fragment07(  5, 3,  4, 16.66*keV ); // Li5
G4StableFermiFragment G4FermiConfiguration::Fragment08(  6, 2,  1,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment09(  6, 3,  3,  0.00*keV );
			                    
G4StableFermiFragment G4FermiConfiguration::Fragment10(  6, 3,  1,  3.56*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment11(  7, 3,  4,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment12(  7, 3,  2,  0.48*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment13(  7, 4,  4,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment14(  7, 4,  2,  0.43*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment15(  8, 3,  5,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment16(  8, 3,  3,  0.98*keV );
G4Be8FermiFragment    G4FermiConfiguration::Fragment17(  8, 4,  1,  0.00*keV ); // Be8
G4StableFermiFragment G4FermiConfiguration::Fragment18(  9, 4,  4,  0.00*keV );
G4B9FermiFragment     G4FermiConfiguration::Fragment19(  9, 5,  4,  0.00*keV ); // B9
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment20( 10, 4,  1,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment21( 10, 4,  5,  3.37*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment22( 10, 4,  8,  5.96*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment23( 10, 4,  1,  6.18*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment24( 10, 4,  5,  6.26*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment25( 10, 5,  7,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment26( 10, 5,  3,  0.72*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment27( 10, 5,  1,  1.74*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment28( 10, 5,  3,  2.15*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment29( 10, 5,  5,  3.59*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment30( 10, 6,  3,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment31( 10, 6,  5,  3.35*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment32( 11, 5,  4,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment33( 11, 5,  2,  2.13*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment34( 11, 5,  6,  4.44*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment35( 11, 5,  4,  5.02*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment36( 11, 5, 10,  6.76*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment37( 11, 5,  6,  7.29*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment38( 11, 5,  4,  7.98*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment39( 11, 5,  6,  8.56*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment40( 11, 6,  4,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment41( 11, 6,  2,  2.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment42( 11, 6,  6,  4.32*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment43( 11, 6,  4,  4.80*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment44( 11, 6,  2,  6.34*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment45( 11, 6,  8,  6.48*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment46( 11, 6,  6,  6.90*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment47( 11, 6,  4,  7.50*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment48( 11, 6,  4,  8.10*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment49( 11, 6,  6,  8.42*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment50( 11, 6,  8,  8.66*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment51( 12, 5,  3,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment52( 12, 5,  5,  0.95*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment53( 12, 5,  5,  1.67*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment54( 12, 5,  4,  2.65*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment55( 12, 6,  1,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment56( 12, 6,  5,  4.44*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment57( 13, 6,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment58( 13, 6,  2,  3.09*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment59( 13, 6,  4,  3.68*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment60( 13, 6,  6,  3.85*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment61( 13, 7,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment62( 14, 6,  1,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment63( 14, 6,  3,  6.09*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment64( 14, 6,  8,  6.69*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment65( 14, 6,  6,  6.96*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment66( 14, 6,  5,  7.34*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment67( 14, 7,  3,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment68( 14, 7,  1,  2.31*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment69( 14, 7,  3,  3.95*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment70( 14, 7,  1,  4.92*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment71( 14, 7,  5,  5.11*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment72( 14, 7,  3,  5.69*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment73( 14, 7,  7,  5.83*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment74( 14, 7,  3,  6.20*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment75( 14, 7,  7,  6.44*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment76( 14, 7,  5,  7.03*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment77( 15, 7,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment78( 15, 7,  8,  5.28*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment79( 15, 7,  4,  6.32*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment80( 15, 7, 10,  7.22*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment81( 15, 7,  8,  7.57*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment82( 15, 7,  2,  8.31*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment83( 15, 7,  4,  8.57*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment84( 15, 7, 14,  9.15*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment85( 15, 7, 14,  9.79*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment86( 15, 7,  8, 10.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment87( 15, 8,  2,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment88( 15, 8,  8,  5.22*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment89( 15, 8,  4,  6.18*keV );
			                     
G4StableFermiFragment G4FermiConfiguration::Fragment90( 15, 8, 10,  6.83*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment91( 15, 8,  8,  7.28*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment92( 16, 7,  5,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment93( 16, 7,  1,  0.12*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment94( 16, 7,  7,  0.30*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment95( 16, 7,  3,  0.40*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment96( 16, 8,  1,  0.00*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment97( 16, 8,  8,  6.10*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment98( 16, 8,  5,  6.92*keV );
G4StableFermiFragment G4FermiConfiguration::Fragment99( 16, 8,  3,  7.12*keV );


                   
G4VFermiFragment * G4FermiConfiguration::theListOfFragments[NumberOfFragments] = {
  &G4FermiConfiguration::Fragment00, 
  &G4FermiConfiguration::Fragment01, 
  &G4FermiConfiguration::Fragment02, 
  &G4FermiConfiguration::Fragment03, 
  &G4FermiConfiguration::Fragment04, 
  &G4FermiConfiguration::Fragment05, 
  &G4FermiConfiguration::Fragment06, 
  &G4FermiConfiguration::Fragment07, 
  &G4FermiConfiguration::Fragment08, 
  &G4FermiConfiguration::Fragment09, 

  &G4FermiConfiguration::Fragment10, 
  &G4FermiConfiguration::Fragment11, 
  &G4FermiConfiguration::Fragment12, 
  &G4FermiConfiguration::Fragment13, 
  &G4FermiConfiguration::Fragment14, 
  &G4FermiConfiguration::Fragment15, 
  &G4FermiConfiguration::Fragment16, 
  &G4FermiConfiguration::Fragment17, 
  &G4FermiConfiguration::Fragment18, 
  &G4FermiConfiguration::Fragment19, 

  &G4FermiConfiguration::Fragment20, 
  &G4FermiConfiguration::Fragment21, 
  &G4FermiConfiguration::Fragment22, 
  &G4FermiConfiguration::Fragment23, 
  &G4FermiConfiguration::Fragment24, 
  &G4FermiConfiguration::Fragment25, 
  &G4FermiConfiguration::Fragment26, 
  &G4FermiConfiguration::Fragment27, 
  &G4FermiConfiguration::Fragment28, 
  &G4FermiConfiguration::Fragment29,

  &G4FermiConfiguration::Fragment30, 
  &G4FermiConfiguration::Fragment31, 
  &G4FermiConfiguration::Fragment32, 
  &G4FermiConfiguration::Fragment33, 
  &G4FermiConfiguration::Fragment34,
  &G4FermiConfiguration::Fragment35, 
  &G4FermiConfiguration::Fragment36, 
  &G4FermiConfiguration::Fragment37, 
  &G4FermiConfiguration::Fragment38, 
  &G4FermiConfiguration::Fragment39,

  &G4FermiConfiguration::Fragment40, 
  &G4FermiConfiguration::Fragment41, 
  &G4FermiConfiguration::Fragment42, 
  &G4FermiConfiguration::Fragment43, 
  &G4FermiConfiguration::Fragment44,
  &G4FermiConfiguration::Fragment45, 
  &G4FermiConfiguration::Fragment46, 
  &G4FermiConfiguration::Fragment47, 
  &G4FermiConfiguration::Fragment48, 
  &G4FermiConfiguration::Fragment49,

  &G4FermiConfiguration::Fragment50, 
  &G4FermiConfiguration::Fragment51, 
  &G4FermiConfiguration::Fragment52, 
  &G4FermiConfiguration::Fragment53, 
  &G4FermiConfiguration::Fragment54,
  &G4FermiConfiguration::Fragment55, 
  &G4FermiConfiguration::Fragment56, 
  &G4FermiConfiguration::Fragment57, 
  &G4FermiConfiguration::Fragment58, 
  &G4FermiConfiguration::Fragment59,

  &G4FermiConfiguration::Fragment60, 
  &G4FermiConfiguration::Fragment61, 
  &G4FermiConfiguration::Fragment62, 
  &G4FermiConfiguration::Fragment63, 
  &G4FermiConfiguration::Fragment64,
  &G4FermiConfiguration::Fragment65, 
  &G4FermiConfiguration::Fragment66, 
  &G4FermiConfiguration::Fragment67, 
  &G4FermiConfiguration::Fragment68, 
  &G4FermiConfiguration::Fragment69,

  &G4FermiConfiguration::Fragment70, 
  &G4FermiConfiguration::Fragment71, 
  &G4FermiConfiguration::Fragment72, 
  &G4FermiConfiguration::Fragment73, 
  &G4FermiConfiguration::Fragment74,
  &G4FermiConfiguration::Fragment75, 
  &G4FermiConfiguration::Fragment76, 
  &G4FermiConfiguration::Fragment77, 
  &G4FermiConfiguration::Fragment78, 
  &G4FermiConfiguration::Fragment79,

  &G4FermiConfiguration::Fragment80, 
  &G4FermiConfiguration::Fragment81, 
  &G4FermiConfiguration::Fragment82, 
  &G4FermiConfiguration::Fragment83, 
  &G4FermiConfiguration::Fragment84,
  &G4FermiConfiguration::Fragment85, 
  &G4FermiConfiguration::Fragment86, 
  &G4FermiConfiguration::Fragment87, 
  &G4FermiConfiguration::Fragment88, 
  &G4FermiConfiguration::Fragment89,

  &G4FermiConfiguration::Fragment90, 
  &G4FermiConfiguration::Fragment91, 
  &G4FermiConfiguration::Fragment92, 
  &G4FermiConfiguration::Fragment93, 
  &G4FermiConfiguration::Fragment94,
  &G4FermiConfiguration::Fragment95, 
  &G4FermiConfiguration::Fragment96, 
  &G4FermiConfiguration::Fragment97, 
  &G4FermiConfiguration::Fragment98, 
  &G4FermiConfiguration::Fragment99
};





G4FermiConfiguration::G4FermiConfiguration()
{
}

G4FermiConfiguration::G4FermiConfiguration(const G4FermiConfiguration &right)
{
  Index = right.Index;
}


G4FermiConfiguration::~G4FermiConfiguration()
{
}


const G4FermiConfiguration & G4FermiConfiguration::operator=(const G4FermiConfiguration &right)
{
  Index = right.Index;
  return *this;
}


G4bool G4FermiConfiguration::operator==(const G4FermiConfiguration &right) const
{
  if (Index.entries() == right.Index.entries()) {
    for (G4int i = 0; i < Index.entries(); i++) {
      if (Index(i) != right.Index(i)) return false;
    }
    return true;
  }
  else return false;
}

G4bool G4FermiConfiguration::operator!=(const G4FermiConfiguration &right) const
{
  return !(*this == right);
}


void G4FermiConfiguration::Initialize(const G4int max)
{
  Index.clear();
  for (G4int i = 0;  i < max; i++) Index.insert(1);
}


G4bool G4FermiConfiguration::SplitNucleus(const G4int A, const G4int Z)
{
  // Splits nucleus (A,Z) into K fragments
  // Returns TRUE if splitting is succesful and FALSE in other case
  
  G4int K = Index.entries();


  G4int L = 0;
  G4int SumA = 0, SumZ = 0;
  for (;;) {
    L++;
    if (L < K) {
      Index[L-1]++;
      if (Index[L-1] > Index[L]) {
	Index[L-1] = 1;
	continue;
      } else {
	SumA = 0;
	for (G4int i = 1; i <= K; i++) SumA += theListOfFragments[Index[i-1]-1]->GetA();
	if (SumA > A) {
	  Index[L-1] = 1;
	  continue;
	} else if (SumA < A) {
	  L = 0;
	  continue;
	} else {
	  SumZ = 0;
	  for (G4int i = 1; i <= K; i++) SumZ += theListOfFragments[Index[i-1]-1]->GetZ();
	  if (SumZ != Z) {
	    L = 0;
	    continue;
	  } else {
	    return true;
	  }
	}
      }
    } else {
      Index[L-1]++;
      if (Index[L-1] > 100) {
	return false;
      } else {
	SumA = 0;
	for (G4int i = 1; i <= K; i++) SumA += theListOfFragments[Index[i-1]-1]->GetA();
	if (SumA < A) {
	  L = 0;
	  continue;
	} else if (SumA == A) {
	  SumZ = 0;
	  for (G4int i = 1; i <= K; i++) SumZ += theListOfFragments[Index[i-1]-1]->GetZ();
	  if (SumZ != Z) {
	    L = 0;
	    continue;
	  } else {
	    return true;
	  }
	} else {
	  return false;
	}
      }
    }
  }
}


G4double G4FermiConfiguration::CoulombBarrier(void)
{
	//  Calculates Coulomb Barrier (MeV) for given channel with K fragments.
	const G4double Coef = (3./5.)*1.44*MeV*fermi* pow(1./(1.+Kappa), 1./3.)/r0;

	G4double SumA = 0, SumZ = 0;
	G4double CoulombEnergy = 0.;
	for (G4int i = 0; i < Index.entries(); i++) {
		G4double z = theListOfFragments[Index[i]-1]->GetZ();
		G4double a = theListOfFragments[Index[i]-1]->GetA();
		CoulombEnergy += (z*z) / pow(a, 1./3.);
		SumA += a;
		SumZ += z;
	}
	CoulombEnergy -= SumZ*SumZ/pow(SumA, 1./3.);
	return -Coef * CoulombEnergy;
}




G4double G4FermiConfiguration::DecayProbability(const G4int A, const G4double TotalE)
  // Decay probability  for a given channel with K fragments
{
  // A: Atomic Weight
	// TotalE: Total energy of nucleus (MeV)
  
	G4int K = Index.entries();
	G4int i;
  
	const G4double NucleonMass = 938.0*MeV;
	const G4double DimCoeff = pow(r0*sqrt(NucleonMass)/hbarc,3.0)*Kappa*sqrt(2.0/pi)/3.0;
  
	// Calculation of 1/Gamma(3(n-1)/2)
	G4double InvGammaFunc = 1.0;
	if (K <= 1) InvGammaFunc = 0.0;
	else {
		G4double arg = 3.0*(K-1)/2.0 - 1.0;
		while (arg > 1.1) {
			InvGammaFunc *= arg; 
			arg--;
		}
		
		if ((K-1)%2 == 1) InvGammaFunc *= sqrt(pi)/2.0;
		
		InvGammaFunc = 1.0/InvGammaFunc;
	}
  
  
	G4double DeltaEnergy = TotalE; // MeV
	G4double Weight = 0.;
	G4double ProdAMass = 1.;  
	G4double ProdSpin = 1.;

	for (i = 0; i<K; i++) {
		ProdAMass *= theListOfFragments[Index[i]-1]->GetA();
		// Spin factor S_n
		ProdSpin *= theListOfFragments[Index[i]-1]->GetPolarization();
		DeltaEnergy -= theListOfFragments[Index[i]-1]->GetFragmentMass() + 
							theListOfFragments[Index[i]-1]->GetExcitationEnergy();
	};
  
	// Check that there is enough energy to produce K fragments
	if ((DeltaEnergy -= CoulombBarrier()) <= 0.0) return Weight; // return 0.0
  
	ProdAMass /= A;
	ProdAMass *= sqrt(ProdAMass);
  
	if (K <= 2) {
		Weight = InvGammaFunc*A*DimCoeff*ProdAMass*ProdSpin*sqrt(DeltaEnergy);
		if (Index[0] == Index[1]) Weight *= 0.5; //Permutation factor G_n
	} else {
		G4double Base = A*DimCoeff*DeltaEnergy*sqrt(DeltaEnergy);
		G4double Powered = 1.0;
		G4double PermutationFactor = 1.0;
		for (G4int i = 0; i < K-1; i++) {
			Powered *= Base;
			G4int N = 1;
			for (G4int j = i+1; j<K; j++) if(Index[i] == Index[j]) N++;
			PermutationFactor *= N;
		};
		Weight = Powered*ProdAMass*ProdSpin*InvGammaFunc/(DeltaEnergy*PermutationFactor);
	}
  
	return Weight; 
}


G4FragmentVector * G4FermiConfiguration::GetFragments(const G4Fragment & theNucleus)
{
 
  G4int K = Index.entries();
  
  // Avalaible kinetic energy of system.
  G4double AvalKineticEnergy = theNucleus.GetExcitationEnergy() +
    G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theNucleus.GetZ(),theNucleus.GetA());
  
  G4int i;
  for (i = 0; i < K; i++) 
    AvalKineticEnergy -= theListOfFragments[Index[i]-1]->GetFragmentMass();
  
  
  // Calculate Momenta of K fragments
  G4RWTPtrOrderedVector<G4LorentzVector>* MomentumComponents = 
    FragmentsMomentum(AvalKineticEnergy);
  
  G4FragmentVector * theResult = new G4FragmentVector;
  
  // Go back to the Lab Frame
  for (i = 0; i < K; i++) {
    
    G4LorentzVector FourMomentum(*(MomentumComponents->at(i)));

    
    // Lorentz boost
    FourMomentum.boost(theNucleus.GetMomentum().boostVector());
    
    G4FragmentVector * fragment = theListOfFragments[Index[i]-1]->GetFragment(FourMomentum);
    
    do {
      theResult->insert(fragment->removeFirst());
    } while (fragment->entries() > 0);
    
    delete fragment;
  }
  
  MomentumComponents->clearAndDestroy();
  delete MomentumComponents;
  
  return theResult;
}



G4RWTPtrOrderedVector<G4LorentzVector>* 
G4FermiConfiguration::FragmentsMomentum(G4double KineticEnergy)
{
  // Calculates momentum for K fragments (Kopylov's method of sampling is used)
  // KinetEnergy is the available kinetic energy
  

  G4int K = Index.entries();

  
  G4RWTPtrOrderedVector<G4LorentzVector>* MomentumList = 
    new G4RWTPtrOrderedVector<G4LorentzVector>(K);
  
  G4double AvalaibleMass = 0; 
  for (G4int i=0; i<K; i++) AvalaibleMass += theListOfFragments[Index[i]-1]->GetFragmentMass();
  
  G4double PFragMagCM = 0.0;
  G4double Mass = AvalaibleMass+KineticEnergy;
  G4LorentzVector PFragCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PFragLab(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestLab(0.0,0.0,0.0,Mass);

  for (G4int l = 0; l < K-1; l++) {
    G4int LK = K - l;
    G4double FragMass = theListOfFragments[Index[LK-1]-1]->GetFragmentMass();
    AvalaibleMass -= FragMass;

    if (LK > 2) KineticEnergy *= RNKSI(LK-1); 
    else KineticEnergy = 0.0;

    G4double RestMass = AvalaibleMass + KineticEnergy;

    PFragMagCM = sqrt(
		      abs((Mass*Mass - (FragMass + RestMass)*(FragMass + RestMass))*
			  (Mass*Mass - (FragMass - RestMass)*(FragMass - RestMass)))
		      )/ (2.0*Mass);
    

    // Create a unit vector with a random direction isotropically distributed
    G4ParticleMomentum RandVector(IsotropicVector(PFragMagCM)); 
    
    PFragCM.setVect(RandVector);
    //    PFragCM.setE((Mass*Mass + FragMass*FragMass - RestMass*RestMass)/(2.0*Mass));
    PFragCM.setE(sqrt(RandVector.mag2()+FragMass*FragMass));

    PRestCM.setVect(-RandVector);
    //    PRestCM.setE((Mass*Mass + RestMass*RestMass - FragMass*FragMass)/(2.0*Mass));
    PRestCM.setE(sqrt(RandVector.mag2()+RestMass*RestMass));


    G4ThreeVector BoostV = PRestLab.boostVector();

    PFragLab = PFragCM;
    PFragLab.boost(BoostV);
    PRestLab = PRestCM;
    PRestLab.boost(BoostV);
    
    MomentumList->prepend(new G4LorentzVector(PFragLab));

    Mass = RestMass;
  }
  
  MomentumList->prepend(new G4LorentzVector(PRestLab));
  return MomentumList;
}


G4double G4FermiConfiguration::RNKSI(const G4int K)
{
	G4double csim = (3.0*K-5.0)/(3.0*K-4.0);
	G4double pex = (3.0*K-5.0)/2.0;
	G4double fcsim = sqrt(1.0-csim)*pow(csim,pex);

	G4double csi = 0.0;
	G4double fcsi= 0.0;
	G4double rf = 0.0;
	do {
		csi = G4UniformRand();
		fcsi = sqrt(1.0-csi)*pow(csi,pex);
		rf = fcsim*G4UniformRand();
	} while (rf > fcsi);
	return csi;
}
    
G4ParticleMomentum G4FermiConfiguration::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ParticleMomentum Vector(Magnitude*cos(Phi)*SinTheta,
			    Magnitude*sin(Phi)*SinTheta,
			    Magnitude*CosTheta);
  return Vector;
}
