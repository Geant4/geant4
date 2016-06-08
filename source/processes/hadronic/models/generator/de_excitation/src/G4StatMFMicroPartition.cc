#include "G4StatMFMicroPartition.hh"



// Copy constructor
G4StatMFMicroPartition::G4StatMFMicroPartition(const G4StatMFMicroPartition & right)
{
	G4Exception("G4StatMFMicroPartition::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMicroPartition & G4StatMFMicroPartition::
operator=(const G4StatMFMicroPartition & right)
{
	G4Exception("G4StatMFMicroPartition::operator= meant to not be accessable");
	return *this;
}


G4bool G4StatMFMicroPartition::operator==(const G4StatMFMicroPartition & right) const
{
//	G4Exception("G4StatMFMicroPartition::operator== meant to not be accessable");
	return false;
}
 

G4bool G4StatMFMicroPartition::operator!=(const G4StatMFMicroPartition & right) const
{
//	G4Exception("G4StatMFMicroPartition::operator!= meant to not be accessable");
	return true;
}



void G4StatMFMicroPartition::CoulombFreeEnergy(const G4double anA)
{
	// This Z independent factor in the Coulomb free energy 
	G4double  CoulombConstFactor = 1.0/
			pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0);
	
	CoulombConstFactor = elm_coupling*(3./5.)*(1. - CoulombConstFactor)/
										G4StatMFParameters::Getr0();

	// We use the aproximation Z_f ~ Z/A * A_f
										
	if (anA == 0 || anA == 1) {
			_theCoulombFreeEnergy.insert(CoulombConstFactor*(theZ/theA)*(theZ/theA));
	} else if (anA == 2 || anA == 3 || anA == 4) {
		// Z/A ~ 1/2
		_theCoulombFreeEnergy.insert(CoulombConstFactor*0.5*
												pow(anA,5./3.));
	} else { // anA > 4
		_theCoulombFreeEnergy.insert(CoulombConstFactor*(theZ/theA)*(theZ/theA)*
												pow(anA,5./3.));	
												
	}
}


G4double G4StatMFMicroPartition::GetCoulombEnergy(void)
{
	G4double  CoulombFactor = 1.0/
			pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0);	
			
	G4double CoulombEnergy = elm_coupling*(3./5.)*theZ*theZ*CoulombFactor/
							  		(G4StatMFParameters::Getr0()*pow(theA,1./3.));
									
									
	for (G4int i = 0; i < _thePartition.entries(); i++) 
			CoulombEnergy += _theCoulombFreeEnergy(i) - elm_coupling*(3./5.)*
									(theZ/theA)*(theZ/theA)*pow(_thePartition(i),5./3.)/
									G4StatMFParameters::Getr0();
									
			
	return CoulombEnergy;
}

G4double G4StatMFMicroPartition::GetPartitionEnergy(const G4double T)
{
	G4double  CoulombFactor = 1.0/
			pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0);	
			
	G4double PartitionEnergy = 0.0;
	
	
	// We use the aprox that Z_f ~ Z/A * A_f
	for (G4int i = 0; i < _thePartition.entries(); i++) {
		if (_thePartition(i) == 0 || _thePartition(i) == 1) {
		
			PartitionEnergy += _theCoulombFreeEnergy(i);
			
		} else if (_thePartition(i) == 2) {
			
			PartitionEnergy +=	
				-2.796 // Binding Energy of deuteron ??????
				+ _theCoulombFreeEnergy(i);
				
		} else if (_thePartition(i) == 3) {
		
			PartitionEnergy +=	
				-9.224 // Binding Energy of trtion/He3 ??????
				+ _theCoulombFreeEnergy(i);
				
		} else if (_thePartition(i) == 4) {
		
			PartitionEnergy +=
				-30.11 // Binding Energy of ALPHA ??????
				+ _theCoulombFreeEnergy(i) 
				+ 4.*T*T/InvLevelDensity(4.);
																
		} else {
												
			PartitionEnergy +=
					//Volume term						
					(- G4StatMFParameters::GetE0() + 
					T*T/InvLevelDensity(_thePartition(i)))
					*_thePartition(i) + 
					 
					// Symmetry term
					G4StatMFParameters::GetGamma0()*
					(1.0-2.0*theZ/theA)*(1.0-2.0*theZ/theA)*_thePartition(i) +  
					
					// Surface term
					(G4StatMFParameters::Beta(T) - T*G4StatMFParameters::DBetaDT(T))*
					pow(_thePartition(i),2./3.) +
					
					// Coulomb term 
	     			_theCoulombFreeEnergy(i);
						
		}
	}
	
	PartitionEnergy += elm_coupling*(3./5.)*theZ*theZ*CoulombFactor/
							 (G4StatMFParameters::Getr0()*pow(theA,1./3.))
							 + (3./2.)*T*(_thePartition.entries()-1);
							 
	return PartitionEnergy;
}


G4double G4StatMFMicroPartition::CalcPartitionTemperature(const G4double U,
								const G4double FreeInternalE0)
{
	G4double PartitionEnergy = GetPartitionEnergy(0.0);
	
	// If this happens, T = 0 MeV, which means that probability for this
	// partition will be 0
	if (abs(U + FreeInternalE0 - PartitionEnergy) < 0.003) return -1.0;
	
	// Calculate temperature by midpoint method
	
	// Bracketing the solution
	G4double Ta = 0.001;
	G4double Tb = max(sqrt(8.0*U/theA),0.0012*MeV);
	G4double Tmid = 0.0;
	
	G4double Da = (U + FreeInternalE0 - GetPartitionEnergy(Ta))/U;
	G4double Db = (U + FreeInternalE0 - GetPartitionEnergy(Tb))/U;
	
	while (Da*Db > 0.0) {
		Tb += 0.5*Tb; 	
		Db = (U + FreeInternalE0 - GetPartitionEnergy(Tb))/U;
	}
	
	G4double eps = 1.0e-14*abs(Ta-Tb);

	for (G4int i = 0; i < 1000; i++) {
		Tmid = (Ta+Tb)/2.0;
		if (abs(Ta-Tb) <= eps) return Tmid;
		G4double Dmid = (U + FreeInternalE0 - GetPartitionEnergy(Tmid))/U;
		if (abs(Dmid) < 0.003) return Tmid;
		if (Da*Dmid < 0.0) {
      	Tb = Tmid;
      	Db = Dmid;
    	} else {
      	Ta = Tmid;
      	Da = Dmid;
    	} 
	}
	// if we arrive here the temperature could not be calculated
	G4cerr << "G4StatMFMicroPartition::CalcPartitionTemperature: I can't calculate the temperature"  
			<< endl;
	// and set probability to 0 returning T < 0
	return -1.0;
	
}


G4double G4StatMFMicroPartition::CalcPartitionProbability(const G4double U,
															const G4double FreeInternalE0,
															const G4double SCompound)
{	
	G4double T = CalcPartitionTemperature(U,FreeInternalE0);
	if ( T <= 0.0) return _Probability = 0.0;
	_Temperature = T;
	

	// Factorial of fragment multiplicity
	G4double Fact = 1.0;
	G4int i;
	for (i = 0; i < _thePartition.entries() - 1; i++) {
		G4double f = 1.0;
		for (G4int ii = i+1; i< _thePartition.entries(); i++) 
								if (_thePartition(i) == _thePartition(ii)) f++;
		Fact *= f;
	}
	
	G4double ProbDegeneracy = 1.0;
	G4double ProbA32 = 1.0;	
	
	for (i = 0; i < _thePartition.entries(); i++) {
		ProbDegeneracy *= GetDegeneracyFactor(G4int(_thePartition(i)));
   	ProbA32 *= G4double(_thePartition(i))*sqrt(G4double(_thePartition(i)));
	}
	
	// Compute entropy
	G4double PartitionEntropy = 0.0;
	for (i = 0; i < _thePartition.entries(); i++) {
		// interaction entropy for alpha
		if (_thePartition(i) == 4)  PartitionEntropy += 
				2.0*T*_thePartition(i)/InvLevelDensity(_thePartition(i));
		// interaction entropy for Af > 4
		else if (_thePartition(i) > 4) PartitionEntropy += 
				2.0*T*_thePartition(i)/InvLevelDensity(_thePartition(i)) -
				G4StatMFParameters::DBetaDT(T)*pow(_thePartition(i),2.0/3.0);
	}
	
	// Thermal Wave Lenght = sqrt(2 pi hbar^2 / nucleon_mass T)
	G4double ThermalWaveLenght3 = 16.15*fermi/sqrt(T);
	ThermalWaveLenght3 = ThermalWaveLenght3*ThermalWaveLenght3*ThermalWaveLenght3;
	
	// Translational Entropy
	G4double kappa = (1. + elm_coupling*(pow(_thePartition.entries(),1./3.)-1.0)/(G4StatMFParameters::Getr0()*pow(theA,1./3.)));
	kappa = kappa*kappa*kappa;
	kappa -= 1.;
	G4double V0 = (4./3.)*pi*theA*G4StatMFParameters::Getr0()*G4StatMFParameters::Getr0()*
						G4StatMFParameters::Getr0();
	G4double FreeVolume = kappa*V0;
	G4double TranslationalS = max(0.0, log(ProbA32/Fact) +
			(_thePartition.entries()-1.0)*log(FreeVolume/ThermalWaveLenght3) +
			1.5*(_thePartition.entries()-1.0) - (3./2.)*log(theA));

	PartitionEntropy += log(ProbDegeneracy) + TranslationalS;
	_Entropy = PartitionEntropy;
	
	// And finally compute probability of fragment configuration
	return _Probability = exp(PartitionEntropy-SCompound);
}



G4double G4StatMFMicroPartition::GetDegeneracyFactor(const G4int A)
{
	// Degeneracy factors are statistical factors
	// DegeneracyFactor for nucleon is (2S_n + 1)(2I_n + 1) = 4
	G4double DegFactor = 0;
	if (A > 4) DegFactor = 1.0;
	else if (A == 1) DegFactor = 4.0;     // nucleon
	else if (A == 2) DegFactor = 3.0;     // Deuteron
	else if (A == 3) DegFactor = 2.0+2.0; // Triton + He3
	else if (A == 4) DegFactor = 1.0;     // alpha
	return DegFactor;
}


G4StatMFChannel * G4StatMFMicroPartition::ChooseZ(const G4double A0, const G4double Z0, const G4double MeanT)
	// Gives fragments charges
{
	G4RWTValOrderedVector<G4int> FragmentsZ;

	G4int ZBalance = 0;
	do {
		G4double CC = G4StatMFParameters::GetGamma0()*8.0;
		G4int SumZ = 0;
		for (G4int i = 0; i < _thePartition.length(); i++) {
			G4double ZMean;
			G4double Af = _thePartition(i);
			if (Af > 1.5 && Af < 4.5) ZMean = 0.5*Af;
			else ZMean = Af*Z0/A0;
			G4double ZDispersion = sqrt(Af * MeanT/CC);
			G4int Zf;
			do {
				Zf = G4int(RandGauss::shoot(ZMean,ZDispersion));
			} while (Zf < 0 || Zf > Af);
			FragmentsZ.insert(Zf);
			SumZ += Zf;
		}
		ZBalance = G4int(Z0) - SumZ;
	} while (abs(ZBalance) > 1.1);
	FragmentsZ(0) += ZBalance;
	
	G4StatMFChannel * theChannel = new G4StatMFChannel;
	for (G4int i = 0; i < _thePartition.length(); i++) 
		theChannel->CreateFragment(_thePartition(i),FragmentsZ(i));
	
	
	return theChannel;
}
