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
//
// by V. Lara
// --------------------------------------------------------------------
//
// Modified:
// 25.07.08 I.Pshenichnov (in collaboration with Alexander Botvina and Igor 
//          Mishustin (FIAS, Frankfurt, INR, Moscow and Kurchatov Institute, 
//          Moscow, pshenich@fias.uni-frankfurt.de) fixed infinite loop for 
//          a fagment with Z=A; fixed memory leak

#include "G4StatMFMacroCanonical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"

// constructor
G4StatMFMacroCanonical::G4StatMFMacroCanonical(const G4Fragment & theFragment) 
{

  // Get memory for clusters
  _theClusters.push_back(new G4StatMFMacroNucleon);              // Size 1
  _theClusters.push_back(new G4StatMFMacroBiNucleon);            // Size 2
  _theClusters.push_back(new G4StatMFMacroTriNucleon);           // Size 3
  _theClusters.push_back(new G4StatMFMacroTetraNucleon);         // Size 4
  for (G4int i = 4; i < theFragment.GetA_asInt(); i++)   
    _theClusters.push_back(new G4StatMFMacroMultiNucleon(i+1)); // Size 5 ... A
  
  // Perform class initialization
  Initialize(theFragment);
    
}

// destructor
G4StatMFMacroCanonical::~G4StatMFMacroCanonical() 
{
  // garbage collection
  if (!_theClusters.empty()) 
    {
      std::for_each(_theClusters.begin(),_theClusters.end(),DeleteFragment());
    }
}

// Initialization method
void G4StatMFMacroCanonical::Initialize(const G4Fragment & theFragment) 
{
  
  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  G4double x = 1.0 - 2.0*Z/G4double(A);
  G4Pow* g4calc = G4Pow::GetInstance();
  
  // Free Internal energy at T = 0
  __FreeInternalE0 = A*( -G4StatMFParameters::GetE0() +     // Volume term (for T = 0)
			 G4StatMFParameters::GetGamma0()*x*x) // Symmetry term
    + G4StatMFParameters::GetBeta0()*g4calc->Z23(A) +   // Surface term (for T = 0)
    0.6*elm_coupling*Z*Z/(G4StatMFParameters::Getr0()*     // Coulomb term 
			  g4calc->Z13(A));
  
  CalculateTemperature(theFragment);
  return;
}

void G4StatMFMacroCanonical::CalculateTemperature(const G4Fragment & theFragment)
{
  // Excitation Energy
  G4double U = theFragment.GetExcitationEnergy();
  
  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  
  // Fragment Multiplicity
  G4double FragMult = std::max((1.0+(2.31/MeV)*(U/A - 3.5*MeV))*A/100.0, 2.0);

  // Parameter Kappa
  G4Pow* g4calc = G4Pow::GetInstance();
  _Kappa = (1.0+elm_coupling*(g4calc->A13(FragMult)-1)/
	    (G4StatMFParameters::Getr0()*g4calc->Z13(A)));
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
  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  
  std::vector<G4int> ANumbers(A);
  
  G4double Multiplicity = ChooseA(A,ANumbers);
  
  std::vector<G4int> FragmentsA;
  
  G4int i = 0;  
  for (i = 0; i < A; i++) 
    {
      for (G4int j = 0; j < ANumbers[i]; j++) FragmentsA.push_back(i+1);
    }
  
  // Sort fragments in decreasing order
  G4int im = 0;
  for (G4int j = 0; j < Multiplicity; j++) 
    {
      G4int FragmentsAMax = 0;
      im = j;
      for (i = j; i < Multiplicity; i++) 
	{
	  if (FragmentsA[i] <= FragmentsAMax) { continue; }
	  else 
	    {
	      im = i;
	      FragmentsAMax = FragmentsA[im];
	    }
	}	
      if (im != j) 
	{
	  FragmentsA[im] = FragmentsA[j];
	  FragmentsA[j]  = FragmentsAMax;
	}
    }
  return ChooseZ(Z,FragmentsA);
}

G4double G4StatMFMacroCanonical::ChooseA(G4int A, std::vector<G4int> & ANumbers)
  // Determines fragments multiplicities and compute total fragment multiplicity
{
  G4double multiplicity = 0.0;
  G4int i;
    
  std::vector<G4double> AcumMultiplicity;
  AcumMultiplicity.reserve(A);

  AcumMultiplicity.push_back((*(_theClusters.begin()))->GetMeanMultiplicity());
  for (std::vector<G4VStatMFMacroCluster*>::iterator it = _theClusters.begin()+1;
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
    for (i = 0; i < A; i++) ANumbers[i] = 0;
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
      CheckA = A - SumA;
      
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (CheckA > 0);
    
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  } while (CheckA < 0 || std::abs(__MeanMultiplicity - multiplicity) > std::sqrt(__MeanMultiplicity) + 0.5);
  
  return multiplicity;
}

G4StatMFChannel * G4StatMFMacroCanonical::ChooseZ(G4int & Z, 
						  std::vector<G4int> & FragmentsA)
    // 
{
  G4Pow* g4calc = G4Pow::GetInstance();
  std::vector<G4int> FragmentsZ;
  
  G4int DeltaZ = 0;
  G4double CP =  G4StatMFParameters::GetCoulomb();
  G4int multiplicity = (G4int)FragmentsA.size();
  
  do {
    FragmentsZ.clear();
    G4int SumZ = 0;
    for (G4int i = 0; i < multiplicity; ++i) 
      {
	G4int A = FragmentsA[i];
	if (A <= 1) 
	  {
	    G4double RandNumber = G4UniformRand();
	    if (RandNumber < (*_theClusters.begin())->GetZARatio()) 
	      {
		FragmentsZ.push_back(1);
		SumZ += FragmentsZ[i];
	      } 
	    else FragmentsZ.push_back(0);
	  } 
	else 
	  {
	    G4double RandZ;
	    G4double CC = 8.0*G4StatMFParameters::GetGamma0()
	      + 2*CP*g4calc->Z23(FragmentsA[i]);
	    G4double ZMean;
	    if (FragmentsA[i] > 1 && FragmentsA[i] < 5) { ZMean = 0.5*FragmentsA[i]; }
	    else {
	      ZMean = FragmentsA[i]*(4.0*G4StatMFParameters::GetGamma0()
				     + _ChemPotentialNu)/CC; 
	    }
	    G4double ZDispersion = std::sqrt(FragmentsA[i]*__MeanTemperature/CC);
	    G4int z;
	    do 
	      {
		RandZ = G4RandGauss::shoot(ZMean,ZDispersion);
		z = G4lrint(RandZ+0.5);
		// Loop checking, 05-Aug-2015, Vladimir Ivanchenko
	      } while (z < 0 || z > A);
	    FragmentsZ.push_back(z);
	    SumZ += z;
	  }
      }
    DeltaZ = Z - SumZ;
  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  } while (std::abs(DeltaZ) > 1);
    
  // DeltaZ can be 0, 1 or -1
  G4int idx = 0;
  if (DeltaZ < 0.0)
    {
      while (FragmentsZ[idx] < 1) { ++idx; }
    }
  FragmentsZ[idx] += DeltaZ;
  
  G4StatMFChannel * theChannel = new G4StatMFChannel;
  for (G4int i = multiplicity-1; i >= 0; --i) 
    {
      theChannel->CreateFragment(FragmentsA[i],FragmentsZ[i]);
    }
 
  return theChannel;
}
