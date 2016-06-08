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
// $Id: G4StatMFMacroCanonical.cc,v 1.11 2002/06/06 17:57:37 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// by V. Lara
// --------------------------------------------------------------------

#include "G4StatMFMacroCanonical.hh"


// constructor
G4StatMFMacroCanonical::G4StatMFMacroCanonical(const G4Fragment & theFragment) 
{

    // Get memory for clusters
    _theClusters.push_back(new G4StatMFMacroNucleon);              // Size 1
    _theClusters.push_back(new G4StatMFMacroBiNucleon);            // Size 2
    _theClusters.push_back(new G4StatMFMacroTriNucleon);           // Size 3
    _theClusters.push_back(new G4StatMFMacroTetraNucleon);         // Size 4
    for (G4int i = 4; i < theFragment.GetA(); i++)   
	_theClusters.push_back(new G4StatMFMacroMultiNucleon(i+1)); // Size 5 ... A
		
    // Perform class initialization
    Initialize(theFragment);

}


// destructor
G4StatMFMacroCanonical::~G4StatMFMacroCanonical() 
{
  // garbage collection
  if (!_theClusters.empty()) {
    G4std::for_each(_theClusters.begin(),_theClusters.end(),DeleteFragment());
  }
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
  
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();

  // Free Internal energy at T = 0
  __FreeInternalE0 = A*( -G4StatMFParameters::GetE0() +          // Volume term (for T = 0)
			 G4StatMFParameters::GetGamma0()*        // Symmetry term
			 (1.0-2.0*Z/A)*(1.0-2.0*Z/A) ) +
    G4StatMFParameters::GetBeta0()*pow(A,2.0/3.0) +              // Surface term (for T = 0)
    (3.0/5.0)*elm_coupling*Z*Z/(G4StatMFParameters::Getr0()*     // Coulomb term 
				pow(A,1.0/3.0));
  
  
  
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
  
  G4std::vector<G4double> ANumbers(A);

  G4double Multiplicity = ChooseA(A,ANumbers);


  G4std::vector<G4double> FragmentsA;
  
    G4int i = 0;  
    for (i = 0; i < A; i++) 
      {
	for (G4int j = 0; j < ANumbers[i]; j++) FragmentsA.push_back(i+1);
      }


    // Sort fragments in decreasing order
    G4int im = 0;
    for (G4int j = 0; j < Multiplicity; j++) 
      {
	G4double FragmentsAMax = 0.0;
	im = j;
	for (i = j; i < Multiplicity; i++) 
	  {
	    if (FragmentsA[i] <= FragmentsAMax) continue;
	    else 
	      {
		im = i;
		FragmentsAMax = FragmentsA[im];
	      }
	  }
	
	if (im != j) 
	  {
	    FragmentsA[im] = FragmentsA[j];
	    FragmentsA[j] = FragmentsAMax;
	  }
      }
    
    return ChooseZ(Z,FragmentsA);
}



G4double G4StatMFMacroCanonical::ChooseA(const G4double A, G4std::vector<G4double> & ANumbers)
  // Determines fragments multiplicities and compute total fragment multiplicity
{
  G4double multiplicity = 0.0;
  G4int i;
  
  
  G4std::vector<G4double> AcumMultiplicity;
  AcumMultiplicity.reserve(A);

  AcumMultiplicity.push_back((*(_theClusters.begin()))->GetMeanMultiplicity());
  for (G4std::vector<G4VStatMFMacroCluster*>::iterator it = _theClusters.begin()+1;
       it != _theClusters.end(); ++it)
    {
      AcumMultiplicity.push_back((*it)->GetMeanMultiplicity()+AcumMultiplicity.back());
    }
  
  G4int CheckA;
  do {
    CheckA = -1;
    G4int SumA = 0;
    G4int ThisOne = 0;
    multiplicity = 0.0;
    for (i = 0; i < A; i++) ANumbers[i] = 0.0;
    do {
      G4double RandNumber = G4UniformRand()*__MeanMultiplicity;
      for (i = 0; i < A; i++) {
	if (RandNumber < AcumMultiplicity[i]) {
	  ThisOne = i;
	  break;
	}
      }
      multiplicity++;
      ANumbers[ThisOne] = ANumbers[ThisOne]+1;
      SumA += ThisOne+1;
      CheckA = G4int(A) - SumA;
      
    } while (CheckA > 0);
    
  } while (CheckA < 0 || abs(__MeanMultiplicity - multiplicity) > sqrt(__MeanMultiplicity) + 1./2.);
  
  return multiplicity;
}


G4StatMFChannel * G4StatMFMacroCanonical::ChooseZ(const G4int & Z, 
						  G4std::vector<G4double> & FragmentsA)
    // 
{
    G4std::vector<G4double> FragmentsZ;
	
    G4double DeltaZ = 0.0;
    G4double CP = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
	(1.0 - 1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

    G4int multiplicity = FragmentsA.size();

    do 
    {
      FragmentsZ.clear();
      G4int SumZ = 0;
      for (G4int i = 0; i < multiplicity; i++) 
	{
	  G4double A = FragmentsA[i];
	  if (A <= 1.0) 
	    {
	      G4double RandNumber = G4UniformRand();
	      if (RandNumber < (*_theClusters.begin())->GetZARatio()) 
		{
		  FragmentsZ.push_back(1.0);
		  SumZ += G4int(FragmentsZ[i]);
		} 
	      else FragmentsZ.push_back(0.0);
	    } 
	  else 
	    {
	      G4double RandZ;
	      G4double CC = 8.0*G4StatMFParameters::GetGamma0()+2.0*CP*pow(FragmentsA[i],2./3.);
	      G4double ZMean;
	      if (FragmentsA[i] > 1.5 && FragmentsA[i] < 4.5) ZMean = 0.5*FragmentsA[i];
	      else ZMean = FragmentsA[i]*(4.0*G4StatMFParameters::GetGamma0()+_ChemPotentialNu)/CC;
	      G4double ZDispersion = sqrt(FragmentsA[i]*__MeanTemperature/CC);
	      G4int z;
	      do 
		{
		  RandZ = G4RandGauss::shoot(ZMean,ZDispersion);
		  z = G4int(RandZ+0.5);
		} while (z < 0 || z > A);
	      FragmentsZ.push_back(z);
	      SumZ += z;
	    }
	}
      DeltaZ = Z - SumZ;
    } while (abs(DeltaZ) > 1.1);
    
    // DeltaZ can be 0, 1 or -1
    FragmentsZ[0] += DeltaZ;
	
    G4StatMFChannel * theChannel = new G4StatMFChannel;
    for (G4int i = multiplicity-1; i >= 0; i--) 
      theChannel->CreateFragment(FragmentsA[i],FragmentsZ[i]);
    
    
    return theChannel;
}
