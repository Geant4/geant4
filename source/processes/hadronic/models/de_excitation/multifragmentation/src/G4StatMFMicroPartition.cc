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

#include "G4StatMFMicroPartition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

// Copy constructor
G4StatMFMicroPartition::G4StatMFMicroPartition(const G4StatMFMicroPartition & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMicroPartition::copy_constructor meant to not be accessible");
}

// Operators

G4StatMFMicroPartition & G4StatMFMicroPartition::
operator=(const G4StatMFMicroPartition & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMicroPartition::operator= meant to not be accessible");
  return *this;
}


G4bool G4StatMFMicroPartition::operator==(const G4StatMFMicroPartition & ) const
{
  //throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMicroPartition::operator== meant to not be accessible");
  return false;
}
 

G4bool G4StatMFMicroPartition::operator!=(const G4StatMFMicroPartition & ) const
{
  //throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMicroPartition::operator!= meant to not be accessible");
  return true;
}

void G4StatMFMicroPartition::CoulombFreeEnergy(G4int anA)
{
  // This Z independent factor in the Coulomb free energy 
  G4double  CoulombConstFactor = G4StatMFParameters::GetCoulomb();

  // We use the aproximation Z_f ~ Z/A * A_f

  G4double ZA = G4double(theZ)/G4double(theA);
										
  if (anA == 0 || anA == 1) 
    {
      _theCoulombFreeEnergy.push_back(CoulombConstFactor*ZA*ZA);
    } 
  else if (anA == 2 || anA == 3 || anA == 4) 
    {
      // Z/A ~ 1/2
      _theCoulombFreeEnergy.push_back(CoulombConstFactor*0.5
				      *anA*G4Pow::GetInstance()->Z23(anA));
    } 
  else  // anA > 4
    {
      _theCoulombFreeEnergy.push_back(CoulombConstFactor*ZA*ZA	
				      *anA*G4Pow::GetInstance()->Z23(anA));
    }
}

G4double G4StatMFMicroPartition::GetCoulombEnergy(void)
{
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double  CoulombFactor = 1.0/g4calc->A13(1.0+G4StatMFParameters::GetKappaCoulomb());	
			
  G4double CoulombEnergy = elm_coupling*0.6*theZ*theZ*CoulombFactor/
    (G4StatMFParameters::Getr0()*g4calc->Z13(theA));
	
  G4double ZA = G4double(theZ)/G4double(theA);
  for (unsigned int i = 0; i < _thePartition.size(); i++) 
    CoulombEnergy += _theCoulombFreeEnergy[i] - elm_coupling*0.6*
      ZA*ZA*_thePartition[i]*g4calc->Z23(_thePartition[i])/
      G4StatMFParameters::Getr0();
		
  return CoulombEnergy;
}

G4double G4StatMFMicroPartition::GetPartitionEnergy(G4double T)
{
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double  CoulombFactor = 1.0/g4calc->A13(1.0+G4StatMFParameters::GetKappaCoulomb());	
  
  G4double PartitionEnergy = 0.0;
  
  // We use the aprox that Z_f ~ Z/A * A_f
  for (unsigned int i = 0; i < _thePartition.size(); i++) 
    {
      if (_thePartition[i] == 0 || _thePartition[i] == 1) 
        {	
          PartitionEnergy += _theCoulombFreeEnergy[i];
        }
      else if (_thePartition[i] == 2) 
        {		
          PartitionEnergy +=	
            -2.796 // Binding Energy of deuteron ??????
            + _theCoulombFreeEnergy[i];		
	}
      else if (_thePartition[i] == 3) 
        {	
          PartitionEnergy +=	
            -9.224 // Binding Energy of trtion/He3 ??????
            + _theCoulombFreeEnergy[i];		
	} 
      else if (_thePartition[i] == 4) 
        {	
          PartitionEnergy +=
            -30.11 // Binding Energy of ALPHA ??????
            + _theCoulombFreeEnergy[i] 
            + 4.*T*T/InvLevelDensity(4.);
	} 
      else 
        {
          PartitionEnergy +=
            //Volume term						
            (- G4StatMFParameters::GetE0() + 
             T*T/InvLevelDensity(_thePartition[i]))
            *_thePartition[i] + 
            
            // Symmetry term
            G4StatMFParameters::GetGamma0()*
            (1.0-2.0*theZ/theA)*(1.0-2.0*theZ/theA)*_thePartition[i] +  
            
            // Surface term
            (G4StatMFParameters::Beta(T) - T*G4StatMFParameters::DBetaDT(T))*
            g4calc->Z23(_thePartition[i]) +
            
            // Coulomb term 
            _theCoulombFreeEnergy[i];
	}
    }
	
  PartitionEnergy += elm_coupling*0.6*theZ*theZ*CoulombFactor/
    (G4StatMFParameters::Getr0()*g4calc->Z13(theA))
    + 1.5*T*(_thePartition.size()-1);
  
  return PartitionEnergy;
}

G4double G4StatMFMicroPartition::CalcPartitionTemperature(G4double U,
							  G4double FreeInternalE0)
{
  G4double PartitionEnergy = GetPartitionEnergy(0.0);
  
  // If this happens, T = 0 MeV, which means that probability for this
  // partition will be 0
  if (std::abs(U + FreeInternalE0 - PartitionEnergy) < 0.003) {
    return -1.0;
  }
    
  // Calculate temperature by midpoint method
	
  // Bracketing the solution
  G4double Ta = 0.001;
  G4double Tb = std::max(std::sqrt(8.0*U/theA),0.0012*MeV);
  G4double Tmid = 0.0;
  
  G4double Da = (U + FreeInternalE0 - GetPartitionEnergy(Ta))/U;
  G4double Db = (U + FreeInternalE0 - GetPartitionEnergy(Tb))/U;
  
  if (Da*Db < 0.0) {
    G4bool yes = false;
    for (G4int i = 0; i < 1000; ++i) { 
      Tb += 0.5*Tb; 	
      Db = (U + FreeInternalE0 - GetPartitionEnergy(Tb))/U;
      if (Da*Db >= 0.0) {
	yes = true;
	break;
      }
    }
    if (!yes) { return -1.0; }
  }
  G4double eps = 1.0e-10*std::abs(Ta-Tb);
  
  for (G4int i = 0; i < 1000; ++i) 
    {
      Tmid = (Ta+Tb)/2.0;
      if (std::abs(Ta-Tb) <= eps) { return Tmid; }
      G4double Dmid = (U + FreeInternalE0 - GetPartitionEnergy(Tmid))/U;
      if (std::abs(Dmid) < 0.003) { return Tmid; }
      if (Da*Dmid < 0.0) 
        {
          Tb = Tmid;
          Db = Dmid;
        } 
      else 
        {
          Ta = Tmid;
          Da = Dmid;
        } 
    }
  return -1.0;  
}

G4double G4StatMFMicroPartition::CalcPartitionProbability(G4double U,
							  G4double FreeInternalE0,
							  G4double SCompound)
{	
  G4double T = CalcPartitionTemperature(U,FreeInternalE0);
  if ( T <= 0.0) return _Probability = 0.0;
  _Temperature = T;
  
  G4Pow* g4calc = G4Pow::GetInstance();
  G4int n = (G4int)_thePartition.size();
  
  // Factorial of fragment multiplicity
  G4double Fact = g4calc->factorial(n);
	
  G4double ProbDegeneracy = 1.0;
  G4double ProbA32 = 1.0;	
  G4double PartitionEntropy = 0.0;
  G4double db = G4StatMFParameters::DBetaDT(T);
	
  for (G4int i = 0; i < n; ++i) {
    G4int par = _thePartition[i];
    ProbDegeneracy *= GetDegeneracyFactor(par);
    ProbA32 *= _thePartition[i]*std::sqrt((G4double)par);

    // interaction entropy for alpha
    if (par == 4) {
      PartitionEntropy += 2.0 * T * par/InvLevelDensity(par);
    }
    // interaction entropy for Af > 4
    else if (par > 4) {
      PartitionEntropy += 2.0 * T * par/InvLevelDensity(par) - db * g4calc->Z23(par);
    } 
  }
	
  // Thermal Wave Lenght = std::sqrt(2 pi hbar^2 / nucleon_mass T)
  G4double ThermalWaveLenght3 = 16.15*fermi/std::sqrt(T);
  ThermalWaveLenght3 = ThermalWaveLenght3*ThermalWaveLenght3*ThermalWaveLenght3;
  
  // Translational Entropy
  G4double r0 = G4StatMFParameters::Getr0();
  G4double kappa = 1. + elm_coupling*(g4calc->Z13(n) - 1.0)/(r0*g4calc->Z13(theA));
  G4double V0 = (4./3.)*pi*theA*r0*r0*r0;
  G4double FreeVolume = (kappa*kappa*kappa - 1.0)*V0;
  G4double TranslationalS = G4Log(ProbA32/Fact)
    + (n - 1)*G4Log(FreeVolume/ThermalWaveLenght3)
    + 1.5*(n - 1) - 1.5*g4calc->logZ(theA);
  TranslationalS = std::max(TranslationalS, 0.0);
  
  PartitionEntropy += G4Log(ProbDegeneracy) + TranslationalS;
  _Entropy = PartitionEntropy;
	
  // And finally compute probability of fragment configuration
  G4double exponent = std::min(PartitionEntropy - SCompound, 200.);
  return _Probability = G4Exp(exponent);
}

G4double G4StatMFMicroPartition::GetDegeneracyFactor(G4int A)
{
  // Degeneracy factors are statistical factors
  // DegeneracyFactor for nucleon is (2S_n + 1)(2I_n + 1) = 4
  G4double DegFactor = 0;
  if (A > 4) DegFactor = 1.0;
  else if (A == 1) DegFactor = 4.0;     // nucleon
  else if (A == 2) DegFactor = 3.0;     // Deuteron
  else if (A == 3) DegFactor = 4.0;     // Triton + He3
  else if (A == 4) DegFactor = 1.0;     // alpha
  return DegFactor;
}

G4StatMFChannel * G4StatMFMicroPartition::ChooseZ(G4int A0, G4int Z0, G4double MeanT)
// Gives fragments charges
{
  std::vector<G4int> FragmentsZ;
  G4int n = (G4int)_thePartition.size();
  
  G4int ZBalance = 0;
  do {
    G4double CC = G4StatMFParameters::GetGamma0()*8.0;
    G4int SumZ = 0;
    for (G4int i = 0; i < n; ++i) {
      G4double ZMean;
      G4double Af = _thePartition[i];
      if (Af > 1.5 && Af < 4.5) { ZMean = 0.5*Af; }
      else ZMean = Af*Z0/A0;
      G4double ZDispersion = std::sqrt(Af * MeanT/CC);
      G4int Zf;
      do {
	Zf = static_cast<G4int>(G4RandGauss::shoot(ZMean, ZDispersion));
      } 
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
      while (Zf < 0 || Zf > Af);
      FragmentsZ.push_back(Zf);
      SumZ += Zf;
    }
    ZBalance = Z0 - SumZ;
  } 
  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  while (std::abs(ZBalance) > 1);
  FragmentsZ[0] += ZBalance;
	
  G4StatMFChannel * theChannel = new G4StatMFChannel;
  for (G4int i = 0; i < n; ++i) {
    theChannel->CreateFragment(_thePartition[i], FragmentsZ[i]);
  }

  return theChannel;
}
