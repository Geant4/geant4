// $Id: G4StatMFMicrocanonical.cc,v 1.3 1998/11/12 16:19:51 allison Exp $

#include "G4StatMFMicrocanonical.hh"


// constructor
G4StatMFMicrocanonical::G4StatMFMicrocanonical(const G4Fragment & theFragment) 
{

  // Get memory for channels
  for (G4int i = 0; i < theFragment.GetA(); i++) 
    theChannels.insert(new G4StatMFFragment);

  // Perform class initialization
  Initialize(theFragment);

}


// destructor
G4StatMFMicrocanonical::~G4StatMFMicrocanonical() 
{
  // garbage collection
  theChannels.clearAndDestroy();
  W2.clear();
  W3.clear();
  W4.clear();

  G4int i;
  for (i = 0; i < ANum2.entries(); i++) delete ANum2(i);
  ANum2.clear();
  for (i = 0; i < ANum3.entries(); i++) delete ANum3(i);
  ANum3.clear();
  for (i = 0; i < ANum4.entries(); i++) delete ANum4(i);
  ANum4.clear();
}


// operators definitions
G4StatMFMicrocanonical & 
G4StatMFMicrocanonical::operator=(const G4StatMFMicrocanonical & right) 
{
  G4Exception("G4StatMFMicrocanonical::operator= meant to not be accessable");
  return *this;
}

G4bool G4StatMFMicrocanonical::operator==(const G4StatMFMicrocanonical & right) const 
{
  return false;
}


G4bool G4StatMFMicrocanonical::operator!=(const G4StatMFMicrocanonical & right) const 
{
  return true;
}


// Initialization method

void G4StatMFMicrocanonical::Initialize(const G4Fragment & theFragment) 
{

  G4int i;
  
  // Excitation Energy (in MeV)
  G4double U = theFragment.GetExcitationEnergy()/MeV;
  
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();


  for (i = 0; i < A; i++) {
    
    // Calculate Inverse Density Levels for channels
    // Epsilon0*(1 + 3 (Af - 1)) ????  --> Ask to Nikolai
    // and correct in G4StatMFFragment
    theChannels(i)->SetInvLevelDensity(i);

    // Z/A ratios
    theChannels(i)->SetZARatio(Z/A);  // optimzation ?
    
    // Degeneracy factors
    theChannels(i)->SetDegeneracyFactor(1.0);
    
    // Multiplicities
    theChannels(i)->SetMultiplicity(0.0);
  }
  // Z/A ratios are for Af > 4
  theChannels(1)->SetZARatio(0.5);
  theChannels(2)->SetZARatio(0.5);
  theChannels(3)->SetZARatio(0.5);
  // Que pasa  con theChannels(0).SetZARatio(0.5) ???? --> Ask to Nikolai
  
  // Degeneracy factors are statistical factors
  // DegeneracyFactor for nucleon is (2S_n + 1)(2I_n + 1) = 4
  theChannels(0)->SetDegeneracyFactor(4.0); // nucleon
  theChannels(1)->SetDegeneracyFactor(3.0); //   .
  theChannels(2)->SetDegeneracyFactor(4.0); //   .
  theChannels(3)->SetDegeneracyFactor(1.0); // alpha



  // Configuration temperature
  G4double TConfiguration = sqrt(U/(0.125*A));

  // Free internal energy at Temperature T = 0
  FreeInternalE0 = A*( -G4StatMFParameters::GetE0() +  // Volume term (for T = 0)
		       G4StatMFParameters::GetGamma0()*(1.0-2.0*Z/A)*(1.0-2.0*Z/A) ) + // Symmetry term
    G4StatMFParameters::GetBeta0()*pow(A,2.0/3.0) + // Surface term (for T = 0)
    (3.0/5.0)*1.44*Z*Z/(G4StatMFParameters::Getr0()*pow(A,1.0/3.0)); // Coulomb term 
  
  // Statistical weights
  W = 0.0;
  WW2 = 0.0;
  WW3 = 0.0;
  WW4 = 0.0;

  // Number of configurations for breakups with multiplicities 2, 3 and 4
  M2 = 0;
  M3 = 0;
  M4 = 0;


  // Mean breakup multiplicity
  MeanMultiplicity = 0.0;

  // Mean channel temperature
  MeanTemperature = 0.0;

  // Mean channel entropy
  MeanEntropy = 0.0;

  // Calculate entropy of compound nucleus
  G4double SCompoundNucleus = CalcEntropyOfCompoundNucleus(theFragment,TConfiguration);

  // I have to change this:
  // -----------------------
  
  // Statistical weight of compound nucleus
  WCompoundNucleus = 1.0; // exp(SCompoundNucleus - SCompoundNucleus);
  
  W += WCompoundNucleus;
  MeanMultiplicity += 1.0 * WCompoundNucleus;
  MeanTemperature += TConfiguration * WCompoundNucleus;
  MeanEntropy += SCompoundNucleus * WCompoundNucleus;

  theChannels(G4int(A) - 1)->SetMultiplicity(
					     theChannels(G4int(A) - 1)->GetMultiplicity()+
					     WCompoundNucleus);

  // -----------------------

  // Maximal fragment multiplicity allowed in direct simulation
  G4int MaxMult = 4;
  if (A > 110) MaxMult = 3;

  // Keep fragment atomic numbers
  G4int * FragmentAtomicNumbers = new G4int(G4int(A+0.5));


  // We distribute A nucleons between m fragments mantaining the order
  // FragmentAtomicNumbers[m-1]>FragmentAtomicNumbers[m-2]>...>FragmentAtomicNumbers[0]
  // Our initial distribution is 
  // FragmentAtomicNumbers[m-1]=A, FragmentAtomicNumbers[m-2]=0, ..., FragmentAtomicNumbers[0]=0
  for (G4int m = 2; m <= MaxMult; m++) { 
    FragmentAtomicNumbers[m-1] = A;
    G4int M1 = m - 1;
    for (i = 0; i < M1; i++) FragmentAtomicNumbers[i] = 0;
    // We try to distribute A nucleons between m fragments
    // DistributeNucleonsBetweenFragments return true if it is possible 
    // and false if not
    while (DistributeNucleonsBetweenFragments(m,FragmentAtomicNumbers)) {
      // For allowed distributions of nucleons 
      // we calculate the configuration probability
      G4double ConfigurationProbability = CalcFragmentsConfigProbability(theFragment,m,
									FragmentAtomicNumbers,
									SCompoundNucleus);
      // Es la suma de todas las probabilidades para cada configuracion
      W += ConfigurationProbability;
      
      for (G4int j = 0; j < m; j++) theChannels(j)->SetMultiplicity(theChannels(j)->GetMultiplicity()+
								    ConfigurationProbability);
      MeanMultiplicity += m*ConfigurationProbability;
      MeanTemperature += TConfiguration * ConfigurationProbability;
      
      if (m == 2) {
	WW2 += ConfigurationProbability;
	M2 += 1;
	W2.insert(ConfigurationProbability); // W2(M2-1)
	//	G4StatMF1DVector tmp(2);
	//	RWTValVector<G4int> tmp(2);
	G4int * tmp = new G4int[2];
	//	tmp(0) = FragmentAtomicNumbers[m-1];
	tmp[0] = FragmentAtomicNumbers[m-1];
	//	tmp(1) = FragmentAtomicNumbers[m-2];
	tmp[1] = FragmentAtomicNumbers[m-2];
	ANum2.insert(tmp);
      } else if (m == 3) {
	WW3 += ConfigurationProbability;
	M3 += 1;
	W3.insert(ConfigurationProbability); // W3(M3-1)
	G4int * tmp = new G4int[3];
	tmp[0] = FragmentAtomicNumbers[m-1];
	tmp[1] = FragmentAtomicNumbers[m-2];
	tmp[3] = FragmentAtomicNumbers[m-3];
       	ANum3.insert(tmp);
      } else if (m == 4) {
	WW4 += ConfigurationProbability;
	M4 += 1;
	if (M4 > 10000) continue;
	W4.insert(ConfigurationProbability); // W4(M4-1)
	G4int * tmp = new G4int[4];
	tmp[0] = FragmentAtomicNumbers[m-1];
	tmp[1] = FragmentAtomicNumbers[m-2];
	tmp[2] = FragmentAtomicNumbers[m-3];
	tmp[3] = FragmentAtomicNumbers[m-4];
	ANum4.insert(tmp);
      }
    }
  }

  if (M4 > 10000) M4 = 10000;

  // Normalization of statistical weights
  for (i = 0; i < M2; i++) W2(i) = W2(i)/W;
  for (i = 0; i < M3; i++) W3(i) = W3(i)/W;  
  for (i = 0; i < M4; i++) W4(i) = W4(i)/W;

  WW2 /= W;
  WW3 /= W;
  WW4 /= W;
  
  MeanMultiplicity /= W;
  MeanTemperature /= W;
  MeanEntropy /= W;
  WCompoundNucleus /= W;

  // garbage collection
  delete [] FragmentAtomicNumbers;

}





G4double G4StatMFMicrocanonical::CalcFreeInternalEnergy(const G4Fragment & theFragment, const G4double & T)
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();
  G4double A13 = pow(A,1.0/3.0);
  G4int Indx = G4int(A)-1;
  
  G4double VolumeTerm = (-G4StatMFParameters::GetE0()+T*T/theChannels(Indx)->GetInvLevelDensity())*A;

  G4double SymmetryTerm = G4StatMFParameters::GetGamma0()*(1.0-2.0*theChannels(Indx)->GetZARatio())*
    (1.0-2.0*theChannels(Indx)->GetZARatio())*A;

  G4double SurfaceTerm = (Beta(T)-T*DBetaDT(T))*A13*A13;

  G4double CoulombTerm = (3.0/5.0)*1.44*Z*Z/(G4StatMFParameters::Getr0()*A13);

  return VolumeTerm + SymmetryTerm + SurfaceTerm + CoulombTerm;

//   return (-G4StatMFParameters::GetE0()+T*T/theChannels(Indx)->GetInvLevelDensity() +         // volume
// 	  G4StatMFParameters::GetGamma0()*(1.0-2.0*theChannels(Indx)->GetZARatio())*
// 	  (1.0-2.0*theChannels(Indx)->GetZARatio())*A +                                    // symmetry
//  	  (Beta(T)-T*DBetaDT(T))*A13*A13 +                                                  // surface
// 	  (3.0/5.0)*1.44*Z*Z/(G4StatMFParameters::Getr0()*A13));                            // Coulomb
}



G4double G4StatMFMicrocanonical::CalcEntropyOfCompoundNucleus(const G4Fragment & theFragment,G4double & TConf)
  // Calculates Temperature and Entropy of compound nucleus
{
  const G4double A = theFragment.GetA();
  const G4double Z = theFragment.GetZ();
  const G4double U = theFragment.GetExcitationEnergy()/MeV;
  const G4double A13 = pow(A,1.0/3.0);

  G4double Ta = max(sqrt(U/(0.125*A)),0.0012); 
  G4double Tb = Ta;
  
  G4double ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Ta);
  G4double Da = (U+G4StatMFParameters::GetE0()-ECompoundNucleus)/U;
  G4double Db = 0.0;


  // bracketing the solution
  if (Da == 0.0) {
    TConf = Ta;
    return 2*Ta*A/theChannels(G4int(A)-1)->GetInvLevelDensity() - 
      DBetaDT(Ta)*A13*A13;
  } else if (Da < 0.0) {
    do {
      Tb -= 0.5*Tb;
      ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tb);
      Db = (U+G4StatMFParameters::GetE0()-ECompoundNucleus)/U;
    } while (Db < 0.0);
  } else {
    do {
      Tb += 0.5*Tb;
      ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tb);
      Db = (U+G4StatMFParameters::GetE0()-ECompoundNucleus)/U;
    } while (Db > 0.0);
  }
    

  G4double eps = 1.0e-14 * abs(Tb-Ta);

  for (G4int i = 0; i < 1000; i++) {
    G4double Tc = (Ta+Tb)/2.0;
    if (abs(Ta-Tb) <= eps) {
      TConf = Tc;
      return 2*Tc*A/theChannels(G4int(A)-1)->GetInvLevelDensity() - 
	DBetaDT(Tc)*A13*A13;
    }
    ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tc);
    G4double Dc = (U+G4StatMFParameters::GetE0()-ECompoundNucleus)/U;

    if (Dc == 0.0) {
      TConf = Tc;
      return 2*Tc*A/theChannels(G4int(A)-1)->GetInvLevelDensity() - 
	DBetaDT(Tc)*A13*A13;
    }
    
    if (Da*Dc < 0.0) {
      Tb = Tc;
      Db = Dc;
    } else {
      Ta = Tc;
      Da = Dc;
    } 
  }

  
//   const G4double HT = 0.5;
//   G4double H = 0.0;
//   G4int counter = 0;
//   do {
//     G4int id = 0;
//     G4double ECompoundNucleus = CalcFreeInternalEnergy(theFragment,T);
//     G4double D = (U+G4StatMFParameters::GetE0()-ECompoundNucleus)/U;
//     if (abs(D) < 0.003) {
//       TConf = T;
//       return 2*T*A/theChannels(G4int(A)-1)->GetInvLevelDensity() - 
// 	DBetaDT(T)*A13*A13;
//     }
//     if (D <= 0.0) H = -HT;
//     else H = HT;
//     if (D < 0.0) {
//       do {
// 	T += H/pow(2.0,id);
// 	if (T > 0.001) break;
// 	H = HT;
//       } while (id++ < 30);
//     } else {
//       do {
// 	T += H/pow(2.0,id);
// 	if (T > 0.001) break;
// 	H = HT;
//       } while (id++ < 30);
//     
//     if (id >= 30) {
//       G4cerr << "G4StatMFMicrocanoncal::CalcEntropyOfCompoundNucleus: suspecting nucleus";
//       return 0.0;
//     }
//   } while (counter++ <= 120);
//  G4cerr << "G4StatMFMicrocanoncal::CalcEntropyOfCompoundNucleus: suspecting nucleus";

  G4cerr << "G4StatMFMicrocanoncal::CalcEntropyOfCompoundNucleus: I can't calculate the temperature";

  return 0.0;
}


G4bool G4StatMFMicrocanonical::DistributeNucleonsBetweenFragments(const G4int & k,
								  G4int * ANumbers)
  // Distributes A nucleons between k fragments
  // mantaining the order ANumbers[k-1] > ANumbers[k-2] > ... > ANumbers[0]
  // If it is possible returns true. In other case returns false
{
  G4int l = 1;
  while (l < k) {
    G4int tmp = ANumbers[l-1] + ANumbers[k-1];
    ANumbers[l-1] += 1;
    ANumbers[k-1] -= 1;
    if (ANumbers[l-1] > ANumbers[l] || ANumbers[k-2] > ANumbers[k-1]) {
      ANumbers[l-1] = 1;
      ANumbers[k-1] = tmp - 1;
      l++;
    } else 
      return true;
  }
  return false;
}


G4double G4StatMFMicrocanonical::CalcFragmentsConfigProbability(const G4Fragment & theFragment,
								const G4int & M,
								const G4int * ANumbers,
								const G4double & SCompound)
  // Calculates the probability of a fragment configuration, where M is the multiplicity, 
  // ANumbers keeps the fragments atomic numbers and SCompound is the entropy of the 
  // compound nucleus
{

  G4int i;

  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();
  G4double U = theFragment.GetExcitationEnergy();

  // Free volume avalible to the traslational motion of fragment
  // Vf = \kappa*V0
  // V0 is system volume
  G4double V0 = (4.0/3.0)*pi*A*G4StatMFParameters::Getr0()*
    G4StatMFParameters::Getr0()*G4StatMFParameters::Getr0();
  
  G4double kappa = (1.0 + (1.44/(G4StatMFParameters::Getr0()*pow(A,1.0/3.0)))*
		    (pow(M,1.0/3.0) - 1.0));
  kappa = kappa*kappa*kappa; 
  kappa -= 1.0;
  
  G4double FreeVolume = kappa*V0;


  // Factorial of fragment multiplicity
  G4double Fact = 1.0;

  for (i = 1; i < M; i++) {
    G4double f = 1.0;
    for (G4int ii = i+1; i<= M; i++) if (ANumbers[i-1] == ANumbers[ii-1]) f++;
    Fact *= f;
  }

  // Calculate energies
  G4double ProbDegeneracy = 1.0;
  G4double ProbA32 = 1.0;
  G4double EnergyConfiguration = 0.0;
  G4double EnergyCoulomb = 0.0;

  G4double pkp13 = 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0);
  G4double cp = (3./5.)*1.44*(1. - pkp13)/G4StatMFParameters::Getr0();

  G4double * ECOLA = new G4double[G4int(A+0.5)];
  G4double * EA = new G4double[G4int(A+0.5)];

  for (i = 0; i < M; i++) {
    G4int intAf = ANumbers[i] - 1;

    ProbDegeneracy *= theChannels(intAf)->GetDegeneracyFactor();
    ProbA32 *= ANumbers[i]*sqrt(double(ANumbers[i]));

    if (ANumbers[i] == 0 || ANumbers[i] == 1) {
      ECOLA[intAf] = cp*theChannels(intAf)->GetZARatio()*theChannels(intAf)->GetZARatio();
      EA[intAf] = ECOLA[intAf];
    } else {
      ECOLA[intAf] = cp*theChannels(intAf)->GetZARatio()*theChannels(intAf)->GetZARatio()*
	pow(ANumbers[i],5.0/3.0);
      if      (ANumbers[i] == 2) EA[intAf] = -2.796*+ECOLA[intAf];
      else if (ANumbers[i] == 3) EA[intAf] = -9.224*+ECOLA[intAf];
      else if (ANumbers[i] == 4) EA[intAf] = -30.224*+ECOLA[intAf];
      else EA[intAf] = (-G4StatMFParameters::GetE0() + G4StatMFParameters::GetGamma0()*
			(1.0-2.0*theChannels(intAf)->GetZARatio())*
			(1.0-2.0*theChannels(intAf)->GetZARatio()))*ANumbers[i] +
	     G4StatMFParameters::GetBeta0()*pow(ANumbers[i],2.0/3.0) +
	     ECOLA[intAf];
    }
    
    EnergyConfiguration += EA[intAf];
    EnergyCoulomb += ECOLA[intAf];
  }

  EnergyConfiguration += (3.0/5.0)*1.44*Z*Z*pkp13/
    (G4StatMFParameters::Getr0()*pow(A,1.0/3.0));
  EnergyCoulomb += (3.0/5.0)*1.44*Z*Z*pkp13/
    (G4StatMFParameters::Getr0()*pow(A,1.0/3.0));


  for (i = 0; i < M; i++) 
    EnergyCoulomb += -(3.0/5.0)*1.44*theChannels(ANumbers[i]-1)->GetZARatio()*
      theChannels(ANumbers[i]-1)->GetZARatio()*pow(ANumbers[i],5.0/3.0)/
      G4StatMFParameters::Getr0();

  if (U+FreeInternalE0-EnergyConfiguration < 0.003) return 0.0;


  // Calculate temperature by iteration
  G4double T = 0.0;
  G4double Ta = max(sqrt(U/(A*0.125)),0.0012); // initial value
  G4double Tb = Ta;

  G4double EConfiguration = CalcEnergyConfiguration(A,Z,M,ECOLA,EA,ANumbers,Ta);
  G4double Da = (U+G4StatMFParameters::GetE0()-EnergyConfiguration)/U;
  G4double Db = 0.0;


  // bracketing the solution
  if (Da == 0.0) T = Ta;
  else if (Da < 0.0) {
    do {
      Tb -= 0.5*Tb;
      EConfiguration = CalcEnergyConfiguration(A,Z,M,ECOLA,EA,ANumbers,Tb);
      Db = (U+G4StatMFParameters::GetE0()-EConfiguration)/U;
    } while (Db < 0.0);
  } else {
    do {
      Tb += 0.5*Tb;
      EConfiguration = CalcEnergyConfiguration(A,Z,M,ECOLA,EA,ANumbers,Tb);
      Db = (U+G4StatMFParameters::GetE0()-EConfiguration)/U;
    } while (Db > 0.0);
  }

  G4double eps = 1.0e-14*abs(Ta-Tb);

  for ( i = 0; i < 1000; i++) {
    G4double Tc = (Ta+Tb)/2.0;
    if (abs(Ta-Tb) <= eps) {
      T = Tc;
      break;
    }
    EConfiguration = CalcEnergyConfiguration(A,Z,M,ECOLA,EA,ANumbers,Tc);
    G4double Dc = (U+G4StatMFParameters::GetE0()-EConfiguration)/U;

    if (Dc == 0.0) {
      T = Tc;
      break;
    }
    
    if (Da*Dc < 0.0) {
      Tb = Tc;
      Db = Dc;
    } else {
      Ta = Tc;
      Da = Dc;
    } 
  }

  if (i == 1000)
    G4cerr << "G4StatMFMicrocanoncal::CalcFragmentsCongifProbability: I can't calculate the temperature";


//   // Calculate temperature by iteration
//   G4double T = max(sqrt(U/(A*0.125)),0.0012); // initial value

//   G4double HT = 0.5;
//   G4double H = 0.0;

//   G4int id = 0;
//   G4int counter = 0;
  
//   do {
//     EnergyConfiguration = 0.0;
//     for (i = 0; i < M; i++) {
//       G4int intAf = ANumbers[i] - 1;
//       if (ANumbers[i] == 0 || ANumbers[i] == 1) {
// 	EA[intAf] = ECOLA[intAf];
//       } else {
// 	if (ANumbers[i] == 2) EA[intAf] = -2.796*ECOLA[intAf];
// 	if (ANumbers[i] == 3) EA[intAf] = -9.224*ECOLA[intAf];
//         if (ANumbers[i] == 4) EA[intAf] = -30.11*ECOLA[intAf];
// 	else EA[intAf] =  (-G4StatMFParameters::GetE0() + T*T/theChannels(intAf)->GetInvLevelDensity() +
// 			   G4StatMFParameters::GetGamma0()*
// 			   (1.0-2.0*theChannels(intAf)->GetZARatio())*
// 			   (1.0-2.0*theChannels(intAf)->GetZARatio()))*ANumbers[i] +
// 	       (Beta(T)-T*DBetaDT(T))*pow(ANumbers[i],2.0/3.0) +
// 	       ECOLA[intAf];
//       }
//       EnergyConfiguration += EA[intAf];
//     }
//     EnergyConfiguration += (3.0/5.0)*1.44*Z*Z*pkp13/(G4StatMFParameters::Getr0()*pow(A,1.0/3.0)) +
//       1.5*T*(M-1.0);

//     G4double D = (U+G4StatMFParameters::GetE0()-EnergyConfiguration)/U;
//     if (abs(D) < 0.003) break;

//     if (D <= 0.0) H = -HT;
//     else H = HT;

//     if (D < 0.0) {
//       do {
// 	T += H/pow(2.0,id);
// 	if (T > 0.001) break;
// 	H = HT;
//       } while (id++ < 30);      
//     } else {
//       do {
// 	T += H/pow(2.0,id);
// 	if (T > 0.001) break;
// 	H = HT;
//       } while (id++ < 30);
//     }

//     if (id >= 30) {
//       G4cerr << "G4StatMFMicrocanoncal::CalcFragmentsConfigProbability: suspecting nucleus";
//       return 0.0;
//     }

//   } while (counter++ <= 120);


  // Compute entropy
  G4double SConfiguration = 0.0;
  
  for (i = 0; i < M; i++) {
    // interaction entropy for alpha
    if (ANumbers[i] == 4) SConfiguration += 
			    2.0*T*ANumbers[i]/theChannels(ANumbers[i]-1)->GetInvLevelDensity();
    // interaction entropy for Af > 4
    else if (ANumbers[i] > 4) SConfiguration += 
				2.0*T*ANumbers[i]/theChannels(ANumbers[i]-1)->GetInvLevelDensity() -
				DBetaDT(T)*pow(ANumbers[i],2.0/3.0);
  }


  // Thermal Wave Lenght
  G4double ThermalWaveLenght = 16.15/sqrt(15.0);
  // Translational Entropy
  G4double STranslational = max(log(ProbA32/Fact) +
				(M-1.0)*log(FreeVolume/ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght) +
				1.5*(M-1.0) - log(A*sqrt(A)),0.0);

  SConfiguration += log(ProbDegeneracy) + STranslational;
  
  // Garbage collection
  delete [] ECOLA;
  delete [] EA;

  // And finally compute probability of fragment configuration
  return exp(SConfiguration-SCompound);

}



void G4StatMFMicrocanonical::ChooseAandZ(const G4Fragment & theFragment)
  // Choice of fragment atomic numbers and charges 
{
  G4int FragmentMultiplicity = 0;

  G4double RandNumber = G4UniformRand();

  if (RandNumber < WCompoundNucleus) FragmentMultiplicity = 1;
  else if (RandNumber < (WCompoundNucleus + WW2)) FragmentMultiplicity = 2;
  else if (RandNumber < (WCompoundNucleus + WW2 + WW3)) FragmentMultiplicity = 3;
  else FragmentMultiplicity = 4;

  FragmentsA.clear(); 
  FragmentsZ.clear(); 
  
  G4int borrar = FragmentsA.entries();
  borrar = FragmentsZ.entries();
 

  if (FragmentMultiplicity == 1) {
    FragmentsA.insert(theFragment.GetA());
    FragmentsZ.insert(theFragment.GetZ());
    return;
  } else if (FragmentMultiplicity == 2) {
    G4double Wp = WCompoundNucleus;
    G4int j = 0;
    for (G4int i = 0; i < M2; i++) {
      Wp += W2(i);
      if (Wp > RandNumber) break;
      j++;
    }
    FragmentsA.insert(ANum2(j)[0]);
    FragmentsA.insert(ANum2(j)[1]);
  } else if (FragmentMultiplicity == 3) {
    G4double Wp = WCompoundNucleus + WW2;
    G4int j = 0;
    for (G4int i = 0; i < M3; i++) {
      Wp += W3(i);
      if (Wp > RandNumber) break;
      j++;
    }
    FragmentsA.insert(ANum3(j)[0]);
    FragmentsA.insert(ANum3(j)[1]);
    FragmentsA.insert(ANum3(j)[2]);
  } else if (FragmentMultiplicity == 4) {
      G4double Wp = WCompoundNucleus + WW2 + WW3;
      G4int j = 0;
      for (G4int i = 0; i < M4; i++) {
	Wp += W4(i);
	if (Wp > RandNumber) break;
	j++;
      }
      FragmentsA.insert(ANum4(j)[0]);
      FragmentsA.insert(ANum4(j)[1]);
      FragmentsA.insert(ANum4(j)[2]);
      FragmentsA.insert(ANum4(j)[3]);
  } else {
    G4Exception("G4StatMFMicrocanonical::ChooseAandZ: FragmentMultiplicity value not allowed");
  }

  for (G4int v = 0; v < FragmentMultiplicity; v++) FragmentsZ.insert(0);
  borrar = FragmentsA.entries();
  borrar = FragmentsZ.entries();
  ChooseZ(theFragment,FragmentMultiplicity);

  Multiplicity = FragmentMultiplicity;

  return;
}


void G4StatMFMicrocanonical::ChooseZ(const G4Fragment & theFragment, 
				     const G4int & FragmentMultiplicity)
  // Gives fragments charges
{
  G4int ZBalance = 0;
  do {
    G4double CC = G4StatMFParameters::GetGamma0()*0.8;
    G4int SumZ = 0;
    for (G4int i = 0; i < FragmentMultiplicity; i++) {
      G4double ZMean = FragmentsA(i)*theFragment.GetZ()/theFragment.GetA();
      if (FragmentsA(i) > 1.5 && FragmentsA(i) < 4.5) ZMean = 0.5*FragmentsA(i);
      G4double ZDispersion = sqrt(FragmentsA(i)*MeanTemperature/CC);
      G4int Zf;
      do {
	Zf = G4int(RandGauss::shoot(ZMean,ZDispersion)+0.5);
      } while (Zf < 0 || Zf > FragmentsA(i));
      FragmentsZ(i) = Zf;
      SumZ += Zf;
    }
    ZBalance = theFragment.GetZ() - SumZ;
  } while (abs(ZBalance) > 1.1);
  FragmentsZ(0) += ZBalance;
  return;
}


G4double G4StatMFMicrocanonical::CalcEnergyConfiguration(const G4double A, const G4double Z, const G4int M, 
							 G4double * ECOLA, G4double * EA, const G4int * ANumbers, 
							 const G4double T)
{
  const G4double pkp13 = 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0);
  G4double EnergyConfiguration = 0.0;
  for (G4int i = 0; i < M; i++) {
    G4int intAf = ANumbers[i] - 1;
    if (ANumbers[i] == 0 || ANumbers[i] == 1) {
      EA[intAf] = ECOLA[intAf];
    } else {
      if (ANumbers[i] == 2) EA[intAf] = -2.796*ECOLA[intAf];
      if (ANumbers[i] == 3) EA[intAf] = -9.224*ECOLA[intAf];
      if (ANumbers[i] == 4) EA[intAf] = -30.11*ECOLA[intAf];
      else EA[intAf] =  (-G4StatMFParameters::GetE0() + T*T/theChannels(intAf)->GetInvLevelDensity() +
			 G4StatMFParameters::GetGamma0()*
			 (1.0-2.0*theChannels(intAf)->GetZARatio())*
			 (1.0-2.0*theChannels(intAf)->GetZARatio()))*ANumbers[i] +
	     (Beta(T)-T*DBetaDT(T))*pow(ANumbers[i],2.0/3.0) +
	     ECOLA[intAf];
    }
    EnergyConfiguration += EA[intAf];
  }
  EnergyConfiguration += (3.0/5.0)*1.44*Z*Z*pkp13/(G4StatMFParameters::Getr0()*pow(A,1.0/3.0)) +
    1.5*T*(M-1.0);

  return EnergyConfiguration;
}
