// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StatMFMacroCanonical.cc,v 1.6 2000/08/03 08:47:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// --------------------------------------------------------------------

#include "G4StatMFMacroCanonical.hh"


// constructor
G4StatMFMacroCanonical::G4StatMFMacroCanonical(const G4Fragment & theFragment) 
{

	// Get memory for clusters
	_theClusters.insert(new G4StatMFMacroNucleon);              // Size 1
	_theClusters.insert(new G4StatMFMacroBiNucleon);            // Size 2
	_theClusters.insert(new G4StatMFMacroTriNucleon);           // Size 3
	_theClusters.insert(new G4StatMFMacroTetraNucleon);         // Size 4
	for (G4int i = 4; i < theFragment.GetA(); i++)   
		_theClusters.insert(new G4StatMFMacroMultiNucleon(i+1)); // Size 5 ... A
		
	// Perform class initialization
	Initialize(theFragment);

}


// destructor
G4StatMFMacroCanonical::~G4StatMFMacroCanonical() 
{
	// garbage collection
	_theClusters.clearAndDestroy();
}

// operators definitions
G4StatMFMacroCanonical & 
G4StatMFMacroCanonical::operator=(const G4StatMFMacroCanonical & right) 
{
	G4Exception("G4StatMFMacroCanonical::operator= meant to not be accessable");
	return *this;
}

G4bool G4StatMFMacroCanonical::operator==(const G4StatMFMacroCanonical & right) const 
{
	return false;
}


G4bool G4StatMFMacroCanonical::operator!=(const G4StatMFMacroCanonical & right) const 
{
	return true;
}


// Initialization method


void G4StatMFMacroCanonical::Initialize(const G4Fragment & theFragment) 
{
  
	// Excitation Energy
	G4double U = theFragment.GetExcitationEnergy();
  
	G4double A = theFragment.GetA();
	G4double Z = theFragment.GetZ();

  
	// Free Internal energy at T = 0
	__FreeInternalE0 = A*( -G4StatMFParameters::GetE0() +              // Volume term (for T = 0)
									G4StatMFParameters::GetGamma0()*
									(1.0-2.0*Z/A)*(1.0-2.0*Z/A) ) +            // Symmetry term
    				G4StatMFParameters::GetBeta0()*pow(A,2.0/3.0) +        // Surface term (for T = 0)
    				(3.0/5.0)*elm_coupling*Z*Z/(G4StatMFParameters::Getr0()*
					pow(A,1.0/3.0));                                       // Coulomb term 
  


  CalculateTemperature(theFragment);

  return;
}





void G4StatMFMacroCanonical::CalculateTemperature(const G4Fragment & theFragment)
{
	// Excitation Energy
	G4double U = theFragment.GetExcitationEnergy();
  
	G4double A = theFragment.GetA();
	G4double Z = theFragment.GetZ();

	// Fragment Multiplicity
	G4double FragMult = G4std::max((1.0+(2.31/MeV)*(U/A - 3.5*MeV))*A/100.0, 2.0);


	// Parameter Kappa
	_Kappa = (1.0+elm_coupling*(pow(FragMult,1./3.)-1)/
				(G4StatMFParameters::Getr0()*pow(A,1./3.)));
	_Kappa = _Kappa*_Kappa*_Kappa - 1.0;

	
	G4StatMFMacroTemperature * theTemp = new	
		G4StatMFMacroTemperature(A,Z,U,__FreeInternalE0,_Kappa,&_theClusters);

	__MeanTemperature = theTemp->CalcTemperature();
	_ChemPotentialNu = theTemp->GetChemicalPotentialNu();
	_ChemPotentialMu = theTemp->GetChemicalPotentialMu();
	__MeanMultiplicity = theTemp->GetMeanMultiplicity();
	__MeanEntropy = theTemp->GetEntropy();
	
	delete theTemp;			

	return;
}


// --------------------------------------------------------------------------

G4StatMFChannel * G4StatMFMacroCanonical::ChooseAandZ(const G4Fragment &theFragment)
  // Calculate total fragments multiplicity, fragment atomic numbers and charges
{
	G4double A = theFragment.GetA();
	G4double Z = theFragment.GetZ();
  
	G4RWTValVector<G4double> ANumbers(A);

	G4double Multiplicity = ChooseA(A,ANumbers);


	G4RWTValOrderedVector<G4double> FragmentsA;
  
	G4int i = 0;  
	G4int s = 0;
	for (i = 0; i < A; i++) {
//		if (ANumbers(i) == 0) continue;
		for (G4int j = 0; j < ANumbers(i); j++) FragmentsA.insert(i+1);
	}


// Sort fragments in decreasing order
	G4int im = 0;
	for (G4int j = 0; j < Multiplicity; j++) {
		G4double FragmentsAMax = 0.0;
		im = j;
		for (i = j; i < Multiplicity; i++) {
			if (FragmentsA(i) <= FragmentsAMax) continue;
			else {
				im = i;
				FragmentsAMax = FragmentsA(im);
			}
		}
	
		if (im != j) {
			FragmentsA(im) = FragmentsA(j);
			FragmentsA(j) = FragmentsAMax;
		}
	}

	return ChooseZ(Z,FragmentsA);
}



G4double G4StatMFMacroCanonical::ChooseA(const G4double A, G4RWTValVector<G4double> & ANumbers)
  // Determines fragments multiplicities and compute total fragment multiplicity
{
	G4double multiplicity = 0.0;
	G4int i;
  
  
// G4double GH = 0.0;
//  for (i = 0; i < A; i++) GH += _theClusters(i)->GetMeanMultiplicity();
// It is the same that _MeanMultiplicity  
  
  
//  G4double SqrtGH = sqrt(GH) + 0.5;
  
  
	G4RWTValVector<G4double> AcumMultiplicity(A);

	AcumMultiplicity(0) = _theClusters(0)->GetMeanMultiplicity();
	for (i = 1; i < A; i++) AcumMultiplicity(i) = AcumMultiplicity(i-1) + 
													_theClusters(i)->GetMeanMultiplicity();                 
	G4int CheckA;
	do {
		CheckA = -1;
		G4int SumA = 0;
		G4int ThisOne = 0;
		multiplicity = 0.0;
		for (i = 0; i < A; i++) ANumbers(i) = 0.0;
		do {
			G4double RandNumber = G4UniformRand()*__MeanMultiplicity;
			for (i = 0; i < A; i++) {
				if (RandNumber < AcumMultiplicity(i)) {
					ThisOne = i;
					break;
				}
			}
			multiplicity++;
			ANumbers(ThisOne) = ANumbers(ThisOne)+1;
			SumA += ThisOne+1;
			CheckA = G4int(A) - SumA;

		} while (CheckA > 0);

	} while (CheckA < 0 || abs(__MeanMultiplicity - multiplicity) > sqrt(__MeanMultiplicity) + 1./2.);
  
	return multiplicity;
}


G4StatMFChannel * G4StatMFMacroCanonical::ChooseZ(const G4int & Z, 
								G4RWTValOrderedVector<G4double> & FragmentsA)
  // 
{
	G4RWTValOrderedVector<G4double> FragmentsZ;
	
	G4double DeltaZ = 0.0;
	G4double CP = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
					(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

	G4int multiplicity = FragmentsA.entries();

	do {
		G4int SumZ = 0;
		for (G4int i = 0; i < multiplicity; i++) {
			G4double A = FragmentsA(i);
			if (A <= 1.0) {
				G4double RandNumber = G4UniformRand();
				if (RandNumber < _theClusters(0)->GetZARatio()) {
					FragmentsZ.insert(1.0);
					SumZ += G4int(FragmentsZ(i));
				} else FragmentsZ.insert(0.0);
			} else {
				G4double RandZ;
				G4double CC = 8.0*G4StatMFParameters::GetGamma0()+2.0*CP*pow(FragmentsA(i),2./3.);
				G4double ZMean;
				if (FragmentsA(i) > 1.5 && FragmentsA(i) < 4.5) ZMean = 0.5*FragmentsA(i);
				else ZMean = FragmentsA(i)*(4.0*G4StatMFParameters::GetGamma0()+_ChemPotentialNu)/CC;
				G4double ZDispersion = sqrt(FragmentsA(i)*__MeanTemperature/CC);
				G4int z;
				do {
					RandZ = G4RandGauss::shoot(ZMean,ZDispersion);
					z = G4int(RandZ+0.5);
				} while (z < 0 || z > A);
					FragmentsZ(i) = z;
					SumZ += z;
			}
		}
		DeltaZ = Z - SumZ;
	} while (abs(DeltaZ) > 1.1);
	
	// DeltaZ can be 0, 1 or -1
	FragmentsZ(0) += DeltaZ;
	
	G4StatMFChannel * theChannel = new G4StatMFChannel;
	for (G4int i = multiplicity-1; i >= 0; i--) 
		theChannel->CreateFragment(FragmentsA(i),FragmentsZ(i));
	

	return theChannel;
}
