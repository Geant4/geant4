

#include "G4StatMFMacrocanonical.hh"


// constructor
G4StatMFMacrocanonical::G4StatMFMacrocanonical(const G4Fragment & theFragment) 
{

  // Get memory for channels
  for (G4int i = 0; i < theFragment.GetA(); i++) 
    theChannels.insert(new G4StatMFFragment);

  // Perform class initialization
  Initialize(theFragment);

}


// destructor
G4StatMFMacrocanonical::~G4StatMFMacrocanonical() 
{
  // garbage collection
  theChannels.clearAndDestroy();
}

// operators definitions
G4StatMFMacrocanonical & 
G4StatMFMacrocanonical::operator=(const G4StatMFMacrocanonical & right) 
{
  G4Exception("G4StatMFMacrocanonical::operator= meant to not be accessable");
  return *this;
}

G4bool G4StatMFMacrocanonical::operator==(const G4StatMFMacrocanonical & right) const 
{
  return false;
}


G4bool G4StatMFMacrocanonical::operator!=(const G4StatMFMacrocanonical & right) const 
{
  return true;
}


// Initialization method


void G4StatMFMacrocanonical::Initialize(const G4Fragment & theFragment) 
{
  
  // Excitation Energy (in MeV)
  G4double U = theFragment.GetExcitationEnergy()/MeV;
  
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();


  for (G4int i = 0; i < A; i++) {

    // Calculate Inverse Density Levels for channels
    theChannels(i)->SetInvLevelDensity(i);
    
  }
  
  // Free Internal energy at T = 0
  FreeInternalE0 = A*( -G4StatMFParameters::GetE0() +  // Volume term (for T = 0)
		       G4StatMFParameters::GetGamma0()*(1.0-2.0*Z/A)*(1.0-2.0*Z/A) ) + // Symmetry term
    G4StatMFParameters::GetBeta0()*pow(A,2.0/3.0) + // Surface term (for T = 0)
    (3.0/5.0)*1.44*Z*Z/(G4StatMFParameters::Getr0()*pow(A,1.0/3.0)); // Coulomb term 
  


  CalculateTemperature(theFragment);

  return;
}





void G4StatMFMacrocanonical::CalculateTemperature(const G4Fragment & theFragment)
{
  // Excitation Energy (in MeV)
  G4double U = theFragment.GetExcitationEnergy()/MeV;
  
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();

  // Fragment Multiplicity
  G4double FragMult = min((1.0+2.31*(U/A - 3.5))*A/100.0,
			  2.0);

  // Parameter Kappa
  G4double Kappa = (1.0+1.44*(pow(FragMult,1./3.)-1)/(1.17*pow(A,1./3.)));
  Kappa = Kappa*Kappa*Kappa - 1.0;

  // Temperature
  G4double Ta = max(sqrt(U/(0.125*A)),0.0012);
  G4double Tb = Ta;

  G4double ExcitEnergyPerNucleon, TotalMultiplicity;
  MeanTemperature = Ta;
  FragmentsExcitationEnergyAndEntropy(theFragment,Kappa,ExcitEnergyPerNucleon,TotalMultiplicity);
  G4double Da = ((U/A) - ExcitEnergyPerNucleon)/(U/A);
  G4double Db = 0.0;


  // bracketing the solution
  if (Da == 0.0) {
    MeanTemperature = Ta;
    return;
  } else if (Da < 0.0) {
    do {
      Tb -= 0.5*Tb; MeanTemperature = Tb;
      FragmentsExcitationEnergyAndEntropy(theFragment,Kappa,ExcitEnergyPerNucleon,TotalMultiplicity);
      Db = ((U/A) - ExcitEnergyPerNucleon)/(U/A);
    } while (Db < 0.0);
  } else {
    do {
      Tb += 0.5*Tb; MeanTemperature = Tb;
      FragmentsExcitationEnergyAndEntropy(theFragment,Kappa,ExcitEnergyPerNucleon,TotalMultiplicity);
      Db = ((U/A) - ExcitEnergyPerNucleon)/(U/A);
    } while (Db > 0.0);
  }

  G4double eps = 1.0e-14 * abs(Tb-Ta);

  for (G4int i = 0; i < 1000; i++) {
    G4double Tc = (Ta+Tb)/2.0;
    if (abs(Ta-Tb) <= eps) {
      MeanTemperature = Tc;
      return;
    }
    MeanTemperature = Tc;
    FragmentsExcitationEnergyAndEntropy(theFragment,Kappa,ExcitEnergyPerNucleon,TotalMultiplicity);
    G4double Dc = ((U/A) - ExcitEnergyPerNucleon)/(U/A);

    if (Dc == 0.0) {
      MeanTemperature = Tc;
      return;
    }

    if (Da*Dc < 0.0) {
      Tb = Tc;
      Db = Dc;
    } else {
      Ta = Tc;
      Da = Dc;
    }
  }
  G4cerr << "G4StatMFMacrocanoncal::CalculateTemperature: I can't calculate the temperature";

  return;

//   // Temperature
//   G4double T = max(sqrt(U/(0.125*A)),0.0012);

//   const G4double HT = 0.5;
//   G4double H = 0.0;
//   G4int counter = 0;

//   do {
//     G4int id = 0;
    
//     G4double ExcitEnergyPerNucleon, TotalMultiplicity;

//     FragmentsExcitationEnergyAndEntropy(theFragment,Kappa,
// 					ExcitEnergyPerNucleon,
// 					TotalMultiplicity);

//     G4double D = ((U/A) - ExcitEnergyPerNucleon)/(U/A);

//     if (abs(D) < 0.003) {
//       MeanMultiplicity = TotalMultiplicity;
//       return;
//     }
    
//     if (D <= 0.0) H = -HT;
//     else H = HT;

//     if (D <= 0.0) {
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
//       G4cerr << "G4StatMFMacrocanoncal::CalculateTemperature: suspecting nucleus" << endl;
//       return;
//     }

//   } while (counter++ <= 60);

//   G4cerr << "G4StatMFMicrocanoncal::CalculateTemperature: suspecting nucleus" << endl;
//   return;

}




void G4StatMFMacrocanonical::FragmentsExcitationEnergyAndEntropy(const G4Fragment & theFragment,
								 const G4double Kappa,
								 G4double & ExcitEnergyPerNucleon, 
								 G4double & TotalMultiplicity)
  // Calculates excitation energy per nucleon and summed fragment multiplicity and entropy
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();

  // Model Parameters
  G4double A13 = pow(A,1./3.);
  G4double A23 = A13*A13;
  G4double R0 = G4StatMFParameters::Getr0()*A13;
  G4double R = R0*pow(1.0+G4StatMFParameters::GetKappaCoulomb(), 1./3.);

  // Calculate fragment charges over fragment atomic numbers ratios.
  CalculateZARatio(theFragment,Kappa);
  // Compute fragment energies
  FragmentEnergies(theFragment,Kappa); 
  // Compute summed fragment entropy
  MeanEntropy = TotalFragmentsEntropy(A,Kappa);  // En realidad deberia ser TotalEntropy???
  // Average total fragment energy
  G4double AverTotalFragEnergy = 0.0; 

  G4int i;
  for (i = 0; i < A; i++) AverTotalFragEnergy += 
				  theChannels(i)->GetMultiplicity()*theChannels(i)->GetEnergy();
  // Add Coulomb energy
  AverTotalFragEnergy += 0.6*1.44*Z*Z/R;
  // Excitation energy per nucleon
  ExcitEnergyPerNucleon = (AverTotalFragEnergy - FreeInternalE0)/A;
 
  TotalMultiplicity = 0.0;
  for (i = 0; i< A; i++) TotalMultiplicity += theChannels(i)->GetMultiplicity();
  
  return;
  
}


void G4StatMFMacrocanonical::CalculateZARatio(const G4Fragment & theFragment, const G4double & Kappa)
  // This calculates fragment charges over fragment atomic numbers
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();

  G4double CP = (0.6*1.44/G4StatMFParameters::Getr0())*
    (1.0-1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0));

//   RWTValOrderedVector<G4int> C1;
//   RWTValOrderedVector<G4int> C2;
  RWTValVector<G4int> C1; C1.reshape(G4int(A+0.5)); 
  RWTValVector<G4int> C2; C2.reshape(G4int(A+0.5)); 

  G4int i;
  for (i = 4; i < A; i++) {
    G4double a = i+1.0;
    G4double CC = 8.0*G4StatMFParameters::GetGamma0()+2.0*CP*pow(a,2.0/3.0);
    C1(i) = 4.0*G4StatMFParameters::GetGamma0()/CC;
    C2(i) = 1.0/CC;
  }

  ChemPotentialNu = (Z/A)*(8.0*G4StatMFParameters::GetGamma0()+2.0*CP*pow(A,2./3.)) -
    4.0*G4StatMFParameters::GetGamma0();

  G4int K = 0;
  G4int K1 = 0, K2 = 0;
  G4int id = 0;
  const G4double HM = 1.0;
  do {
    theChannels(1)->SetZARatio(0.5);
    theChannels(2)->SetZARatio(0.5);
    theChannels(3)->SetZARatio(0.5);

    for (i = 4; i < A; i++) 
      { G4int temp = C2(i);
      theChannels(i)->SetZARatio(C1(i) + temp * ChemPotentialNu); }

    // Calculate fragment multiplicities
    CalculateMultiplicities(theFragment,Kappa);
    
    theChannels(0)->SetZARatio(YP/(YP+YN));

    G4double ZTotal = 0.0;
    for (i = 0; i < Z; i++) ZTotal += G4double(i)*theChannels(i)->GetZARatio()*
			      theChannels(i)->GetMultiplicity();
    K++;
    
    G4double D = (Z - ZTotal)/Z;
    if (abs(D) < 0.002) break;
    
    G4double H;
    if (D < 0.0) H = -HM;
    else H = HM;

    if (D <= 0.0) {
      K1 = 1;
      if (K1 == 1 && K2 == 2) id++;
      ChemPotentialNu += pow(H/2.0,id);
    } else {
      K2 = 1;
      if (K1 == 1 && K2 == 2) id++;
      ChemPotentialNu += pow(H/2.0,id);
    }
  } while (K <= 60 && id <= 30);

  if (K > 60 || id > 30) G4cerr << "G4StatMFMacrocanonical::CalculateZARatio: suspecting nucleus" << endl;
  return;
}



void G4StatMFMacrocanonical::CalculateMultiplicities(const G4Fragment & theFragment, const G4double & Kappa)
  // 
{
  G4double A = theFragment.GetA();

  G4double CP = (0.6*1.44/G4StatMFParameters::Getr0())*
    (1.0-1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0));

  // calculation of chemical potential \mu
  ChemPotentialMu = -G4StatMFParameters::GetE0()+
    MeanTemperature*MeanTemperature/theChannels(3)->GetInvLevelDensity() -
    ChemPotentialNu*theChannels(3)->GetZARatio() + 
    G4StatMFParameters::GetGamma0()*(1.0-2.0*theChannels(3)->GetZARatio())*(1.0-2.0*theChannels(3)->GetZARatio()) +
    (2.0/3.0)*Beta(MeanTemperature)/1.71 +
    5.0*CP*theChannels(3)->GetZARatio()*theChannels(3)->GetZARatio()*2.92/3.0 -
    1.5*MeanTemperature/5.0;


  G4int K = 0;
  G4int K1 = 0, K2 = 0;
  G4int id = 0;
  G4double H = 0.0;
  const G4double HM = 1.0;
  do {
    // Calculate mean fragment multiplicities
    MeanFragmentMultiplicities(theFragment, Kappa);

    // found fragment multiplicities should satisfy constraint: \sum_f N(f) A_f = A
    // using this constraint, chemical potential \mu is defined by iterations
    G4double Atot = 0.0;

    for (G4int i = 0; i < A; i++) Atot += i*theChannels(i)->GetMultiplicity();

    K++;

    G4double D = (A-Atot)/A;
    if (abs(D) < 0.001) break; // or break;

    if (D < 0.0) H = -HM;
    else H = HM;

    if ( D<= 0.0 ) {
      K1 = 1;
      if (K1 == 1 && K2 == 1) id++;
      ChemPotentialMu += H/pow(2.0,id);
    } else {
      K2 = 1;
      if (K1 == 1 && K2 == 1) id++;
      ChemPotentialMu += H/pow(2.0,id);
    }
  } while (K <= 60 && id <= 30);

  if (K > 60 || id > 30) G4cerr << "G4StatMFMacrocanonical::CalculateMultiplicities: suspecting nucleus" << endl;
  return;
}


void G4StatMFMacrocanonical::MeanFragmentMultiplicities(const G4Fragment & theFragment, const G4double & Kappa)
  // Calculates fragment multiplicities
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();

  G4double A13 = pow(A,1.0/3.0);
  G4double R0 = G4StatMFParameters::Getr0()*A13;
  G4double V0 = (4.0/3.0)*pi*R0*R0*R0;
  G4double CP = (0.6*1.44/G4StatMFParameters::Getr0())*
    (1.0-1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0));

  G4double ThermalWaveLenght = 16.15/sqrt(MeanTemperature);
  G4double ThermalWaveLenght3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
  
  YN = (2.0*Kappa*V0/ThermalWaveLenght3)*
    exp(ChemPotentialMu/MeanTemperature);

  YP = (2.0*Kappa*V0/ThermalWaveLenght3)*
    exp((ChemPotentialMu+ChemPotentialNu-CP)/MeanTemperature);
  
  Y2 = (3.0*2.0*sqrt(2.0)*Kappa*V0/ThermalWaveLenght3)*
    exp((2.796+2.0*(ChemPotentialMu+ChemPotentialNu*theChannels(1)->GetZARatio())-
	 CP*theChannels(1)->GetZARatio()*theChannels(1)->GetZARatio()*
	 pow(2.0,5.0/3.0))/MeanTemperature);

  Y3 = (4.0*3.0*sqrt(3.0)*Kappa*V0/ThermalWaveLenght3)*
    exp((9.224+3.0*(ChemPotentialMu+ChemPotentialNu*theChannels(2)->GetZARatio())-
	 CP*theChannels(2)->GetZARatio()*theChannels(2)->GetZARatio()*
	 pow(3.0,5.0/3.0))/MeanTemperature);


  Y4 = (8.0*Kappa*V0/ThermalWaveLenght3)*
    exp((30.11+4.0*(ChemPotentialMu+ChemPotentialNu*theChannels(3)->GetZARatio()+
		    MeanTemperature*MeanTemperature/theChannels(3)->GetInvLevelDensity())-
	 CP*theChannels(3)->GetZARatio()*theChannels(3)->GetZARatio()*pow(4.0,5.0/3.0))/MeanTemperature);

  // Nucleon multiplicity
  theChannels(0)->SetMultiplicity(YN+YP);
  // A_f = 2 multiplicity
  theChannels(1)->SetMultiplicity(Y2);
  // A_f = 3 multiplicity
  theChannels(2)->SetMultiplicity(Y3);
  // A_f = 4 multiplicity
  theChannels(3)->SetMultiplicity(Y4);


  for (G4int i = 4; i < A; i++) {
    G4double a = G4double(i) + 1.0;
    G4double a23 = pow(a,2.0/3.0);
    G4double VNE = (ChemPotentialMu+ChemPotentialNu*theChannels(i)->GetZARatio()+
		    G4StatMFParameters::GetE0()+
		    MeanTemperature*MeanTemperature/theChannels(i)->GetInvLevelDensity()-
		    G4StatMFParameters::GetGamma0()*(1.0-2.0*theChannels(i)->GetZARatio())*
		    (1.0-2.0*theChannels(i)->GetZARatio()))*a - 
      Beta(MeanTemperature)*a23 - 
      CP*theChannels(i)->GetZARatio()*theChannels(i)->GetZARatio()*a*a23;
    VNE /= MeanTemperature;

    if (VNE > 30.0) theChannels(i)->SetMultiplicity(999.0);
    else {
      G4double VN = exp(VNE);
      VN *= Kappa*V0*sqrt(a)*a/ThermalWaveLenght3;
      if (VN < 1.0e-30) theChannels(i)->SetMultiplicity(0.0);
      else theChannels(i)->SetMultiplicity(VN);
    }
  }
  return;
}


void G4StatMFMacrocanonical::FragmentEnergies(const G4Fragment & theFragment,const G4double & Kappa)
  // Calculate Fragment energies at actual temperature
{
  G4double A = theFragment.GetA();

  G4double PkP13 = pow(1./(1. + G4StatMFParameters::GetKappaCoulomb()),1./3.);
  // factor needed  for calculate Coulomb energy
  G4double CP = (0.6*1.44/G4StatMFParameters::Getr0())*(1.-PkP13);


  theChannels(0)->SetEnergy(CP*theChannels(0)->GetZARatio() + 1.5*MeanTemperature);
  theChannels(1)->SetEnergy(-2.796 + 
			    CP*theChannels(1)->GetZARatio()*theChannels(1)->GetZARatio()*pow(2.,5./3.) +
			    1.5*MeanTemperature);
  theChannels(2)->SetEnergy(-9.224 + 
			    CP*theChannels(2)->GetZARatio()*theChannels(2)->GetZARatio()*pow(3.,5./3.) +
			    1.5*MeanTemperature);
  theChannels(3)->SetEnergy(-30.11 + 
			    CP*theChannels(3)->GetZARatio()*theChannels(3)->GetZARatio()*pow(4.,5./3.) +
			    1.5*MeanTemperature +
			    4.0*MeanTemperature*MeanTemperature/theChannels(3)->GetInvLevelDensity());

  for (G4int i = 4; i < A; i++) {
    G4double a = i+1.0;
    G4double a23 = pow(a,2./3.);
    
    // Volume and symmetry terms
    G4double EV = a*(MeanTemperature*MeanTemperature/theChannels(i)->GetInvLevelDensity() -
		     G4StatMFParameters::GetE0() +
		     G4StatMFParameters::GetGamma0()*
		     (1.-2.*theChannels(i)->GetZARatio()*theChannels(i)->GetZARatio()));

    // Surface term
    G4double ES = (Beta(MeanTemperature) - MeanTemperature*DBetaDT(MeanTemperature))*a23;

    // Coulomb term
    G4double EC = CP*a23*a*theChannels(i)->GetZARatio()*theChannels(i)->GetZARatio();
    
    // translational term
    G4double ET = 1.5*MeanTemperature;


    // Total Energy
    theChannels(i)->SetEnergy( EV + ES + EC + ET );
    
  }
  return;
}


G4double G4StatMFMacrocanonical::TotalFragmentsEntropy(const G4double & A, const G4double & Kappa)
  // Calculates summed fragments entropy
{
  // Thermal Wave Length at actual temperature
  G4double ThermalWaveLength = 16.15/sqrt(MeanTemperature);
  G4double ThermalWaveLength3 = ThermalWaveLength*ThermalWaveLength*ThermalWaveLength;

  G4double R0 = G4StatMFParameters::Getr0()*pow(A,2./3.);
  G4double V0 = (4.*pi/3.)*R0*R0*R0;

  // Entropy
  G4double S = 0.0;

  if (YN > 0.0) S += YN*(2.5+log(2.*Kappa*V0/(ThermalWaveLength3*YN)));
  if (YP > 0.0) S += YP*(2.5+log(2.*Kappa*V0/(ThermalWaveLength3*YP)));
  if (Y2 > 0.0) S += Y2*(2.5+log(3.*Kappa*V0*2.*sqrt(2.)/(ThermalWaveLength3*Y2)));
  if (Y3 > 0.0) S += Y3*(2.5+log(4.*Kappa*V0*3.*sqrt(3.)/(ThermalWaveLength3*Y3)));
  if (Y4 > 0.0) S += Y4*(2.5+log(8.*Kappa*V0/(ThermalWaveLength3*Y4)) +
			 8.0*MeanTemperature/theChannels(3)->GetInvLevelDensity());

  for (G4int i = 4; i < A; i++) {
    if (theChannels(i)->GetMultiplicity() <= 0.0) continue;
    G4double a = G4double(i)+1.0;
    G4double SV = 2.0*a*MeanTemperature/theChannels(i)->GetInvLevelDensity();
    G4double SS = -DBetaDT(MeanTemperature)*pow(a,2./3.);
    G4double ST = 2.5+log(Kappa*V0*sqrt(a)*a/(ThermalWaveLength3*theChannels(i)->GetMultiplicity()));

    S += (SV+SS+ST)*theChannels(i)->GetMultiplicity();
  }
  
  return S;
}


void G4StatMFMacrocanonical::ChooseAandZ(const G4Fragment &theFragment)
  // Calculate total fragments multiplicity, fragment atomic numbers and charges
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();
  
  RWTValVector<G4double> ANumbers(A);

  Multiplicity = ChooseA(A,ANumbers);

  G4int i;
  for (i = 0; i < Multiplicity; i++) {
    FragmentsA.insert(0.0);
    FragmentsZ.insert(0.0);
  }
  
  G4int s = 0;
  for (i = 0; i < A; i++) {
    if (ANumbers(i) == 0) continue;
    for (G4int j = 0; j < ANumbers(i); j++) FragmentsA(s+j) = i+1;
    s += ANumbers(i) - 1;
  }
  

  G4int im = 0;
  for (G4int j = 0; j < Multiplicity; j++) {
    G4double FragmentsAMax = 0.0;
    for (i = j; i < Multiplicity; i++) {
      if (FragmentsA(i) <= FragmentsAMax) continue;
      else {
	im = i;
	FragmentsAMax = FragmentsA(im);
      }
    }
    FragmentsA(im) = FragmentsA(j);
    FragmentsA(j) = FragmentsAMax;
  }

  ChooseZ(Z,Multiplicity);

  return;

}


G4double G4StatMFMacrocanonical::ChooseA(const G4double A, RWTValVector<G4double> & ANumbers)
  // Determines fragments multiplicities and compute total fragment multiplicity
{
  G4double multiplicity = 0.0;
  G4double GH = 0.0;
  G4int i;
  for (i = 0; i < A; i++) GH += theChannels(i)->GetMultiplicity();
  
  G4double SqrtGH = sqrt(GH) + 0.5;
  RWTValVector<G4double> AcumMultiplicity;

  AcumMultiplicity(0) = theChannels(0)->GetMultiplicity();
  for (i = 1; i < A; i++) AcumMultiplicity(i) = AcumMultiplicity(i-1) + theChannels.at(i)->GetMultiplicity();

  do {
    G4int CheckA = -1;
    G4int SumA = 0;
    G4int ThisOne = 0;
    do {
      if (CheckA < 0) {
	SumA = 0;
	ThisOne = 0;
	for (i = 0; i < A; i++) ANumbers(i) = 0.0;
	multiplicity = 0.0;
      }
      G4double RandNumber = G4UniformRand()*GH;
      for (i = 0; i < A; i++) {
	if (RandNumber < AcumMultiplicity(i)) {
	  ThisOne = i;
	  break;
	}
      }
      multiplicity++;
      ANumbers(ThisOne) = ANumbers(ThisOne)+1;
      SumA += ThisOne+1;
      CheckA = A - SumA;

    } while (CheckA != 0);

  } while (abs(GH - multiplicity) > SqrtGH);
  
  return multiplicity;
}


void G4StatMFMacrocanonical::ChooseZ(const G4int & Z, const G4double multiplicity)
  // 
{
  G4double DeltaZ = 0.0;
  G4double CP = (0.6*1.44/G4StatMFParameters::Getr0())*
    (1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

  do {
    G4int SumZ = 0;
    for (G4int i = 0; i < multiplicity; i++) {
      G4double A = FragmentsA(i);
      if (A <= 1.0) {
	G4double RandNumber = G4UniformRand();
	FragmentsZ.insert(0.0);   // FragmentsZ(i) = 0.0;
	if (RandNumber > (YN/(YN+YP))) {
	  FragmentsZ(i) = FragmentsZ(i) + 1.;
	  SumZ += FragmentsZ(i);
	}
      } else {
	G4double RandZ;
	G4double CC = 8.0*G4StatMFParameters::GetGamma0()+2.0*CP*pow(FragmentsA(i),2./3.);
	G4double ZMean;
	if (FragmentsA(i) > 1.5 && FragmentsA(i) < 4.5) ZMean = 0.5*FragmentsA(i);
	else ZMean = FragmentsZ(i)*(4.*G4StatMFParameters::GetGamma0()+ChemPotentialNu)/CC;
	G4double ZDispersion = sqrt(FragmentsA(i)*MeanTemperature/CC);
	G4int z;
	do {
	  RandZ = RandGauss::shoot(ZMean,ZDispersion);
	  z = G4int(RandZ);
	} while (z < 0 || z > Z);
	FragmentsZ(i) = RandZ;
	SumZ += z;
      }
    }
    DeltaZ = Z - SumZ;
  } while (abs(DeltaZ) > 1.1);
  // DeltaZ can be 0, 1 or -1
  FragmentsZ(0) += DeltaZ;
  return;
}
