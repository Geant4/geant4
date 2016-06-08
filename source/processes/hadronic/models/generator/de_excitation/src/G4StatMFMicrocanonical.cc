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
// $Id: G4StatMFMicroCanonical.cc,v 1.10 2002/06/06 17:57:46 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFMicroCanonical.hh"
#include <numeric>


// Copy constructor
G4StatMFMicroCanonical::G4StatMFMicroCanonical(const G4StatMFMicroCanonical &
					       right)
{
    G4Exception("G4StatMFMicroCanonical::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMicroCanonical & G4StatMFMicroCanonical::
operator=(const G4StatMFMicroCanonical & right)
{
    G4Exception("G4StatMFMicroCanonical::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMicroCanonical::operator==(const G4StatMFMicroCanonical & right) const
{
    G4Exception("G4StatMFMicroCanonical::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMicroCanonical::operator!=(const G4StatMFMicroCanonical & right) const
{
    G4Exception("G4StatMFMicroCanonical::operator!= meant to not be accessable");
    return true;
}



// constructor
G4StatMFMicroCanonical::G4StatMFMicroCanonical(G4Fragment const & theFragment) 
{
    // Perform class initialization
    Initialize(theFragment);

}


// destructor
G4StatMFMicroCanonical::~G4StatMFMicroCanonical() 
{
  // garbage collection
  if (!_ThePartitionManagerVector.empty()) {
    G4std::for_each(_ThePartitionManagerVector.begin(),
		    _ThePartitionManagerVector.end(),
		    DeleteFragment());
  }
}



// Initialization method

void G4StatMFMicroCanonical::Initialize(const G4Fragment & theFragment) 
{
  
  G4std::vector<G4StatMFMicroManager*>::iterator it;

  // Excitation Energy 
  G4double U = theFragment.GetExcitationEnergy();
  
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();
  
  // Configuration temperature
  G4double TConfiguration = sqrt(8.0*U/A);
  
  // Free internal energy at Temperature T = 0
  __FreeInternalE0 = A*( 
			// Volume term (for T = 0)
			-G4StatMFParameters::GetE0() +  
			// Symmetry term
			G4StatMFParameters::GetGamma0()*(1.0-2.0*Z/A)*(1.0-2.0*Z/A) 
			) + 
    // Surface term (for T = 0)
    G4StatMFParameters::GetBeta0()*pow(A,2.0/3.0) + 
    // Coulomb term 
    elm_coupling*(3.0/5.0)*Z*Z/(G4StatMFParameters::Getr0()*pow(A,1.0/3.0));
  
  // Statistical weight
  G4double W = 0.0;
  
  
  // Mean breakup multiplicity
  __MeanMultiplicity = 0.0;
  
  // Mean channel temperature
  __MeanTemperature = 0.0;
  
  // Mean channel entropy
  __MeanEntropy = 0.0;
  
  // Calculate entropy of compound nucleus
  G4double SCompoundNucleus = CalcEntropyOfCompoundNucleus(theFragment,TConfiguration);
  
  // Statistical weight of compound nucleus
  _WCompoundNucleus = 1.0; // exp(SCompoundNucleus - SCompoundNucleus);
  
  W += _WCompoundNucleus;
  
  
  
  // Maximal fragment multiplicity allowed in direct simulation
  G4int MaxMult = G4StatMFMicroCanonical::MaxAllowedMultiplicity;
  if (A > 110) MaxMult -= 1;
  
  
  
  for (G4int m = 2; m <= MaxMult; m++) {
    G4StatMFMicroManager * aMicroManager = 
      new G4StatMFMicroManager(theFragment,m,__FreeInternalE0,SCompoundNucleus);
    _ThePartitionManagerVector.push_back(aMicroManager);
  }
  
  // W is the total probability
  W = G4std::accumulate(_ThePartitionManagerVector.begin(),
			_ThePartitionManagerVector.end(),
			W,SumProbabilities());
  
  // Normalization of statistical weights
  for (it =  _ThePartitionManagerVector.begin(); it !=  _ThePartitionManagerVector.end(); ++it) 
    {
      (*it)->Normalize(W);
    }

  _WCompoundNucleus /= W;
  
  __MeanMultiplicity += 1.0 * _WCompoundNucleus;
  __MeanTemperature += TConfiguration * _WCompoundNucleus;
  __MeanEntropy += SCompoundNucleus * _WCompoundNucleus;
  

  for (it =  _ThePartitionManagerVector.begin(); it !=  _ThePartitionManagerVector.end(); ++it) 
    {
      __MeanMultiplicity += (*it)->GetMeanMultiplicity();
      __MeanTemperature += (*it)->GetMeanTemperature();
      __MeanEntropy += (*it)->GetMeanEntropy();
    }
  
  return;
}





G4double G4StatMFMicroCanonical::CalcFreeInternalEnergy(const G4Fragment & theFragment, 
							const G4double T)
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();
  G4double A13 = pow(A,1.0/3.0);
  
  G4double InvLevelDensityPar = G4StatMFParameters::GetEpsilon0()*(1.0 + 3.0/(A-1.0));
  
  G4double VolumeTerm = (-G4StatMFParameters::GetE0()+T*T/InvLevelDensityPar)*A;
  
  G4double SymmetryTerm = G4StatMFParameters::GetGamma0()*(A - 2.0*Z)*(A - 2.0*Z)/A;
  
  G4double SurfaceTerm = (G4StatMFParameters::Beta(T)-T*G4StatMFParameters::DBetaDT(T))*A13*A13;
  
  G4double CoulombTerm = elm_coupling*(3.0/5.0)*Z*Z/(G4StatMFParameters::Getr0()*A13);
  
  return VolumeTerm + SymmetryTerm + SurfaceTerm + CoulombTerm;
}



G4double G4StatMFMicroCanonical::CalcEntropyOfCompoundNucleus(const G4Fragment & theFragment,G4double & TConf)
  // Calculates Temperature and Entropy of compound nucleus
{
  const G4double A = theFragment.GetA();
  //    const G4double Z = theFragment.GetZ();
  const G4double U = theFragment.GetExcitationEnergy();
  const G4double A13 = pow(A,1.0/3.0);
  
  G4double Ta = G4std::max(sqrt(U/(0.125*A)),0.0012*MeV); 
  G4double Tb = Ta;
  
  G4double ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Ta);
  G4double Da = (U+__FreeInternalE0-ECompoundNucleus)/U;
  G4double Db = 0.0;
  
  
  G4double InvLevelDensity = CalcInvLevelDensity(A);
  
  
  // bracketing the solution
  if (Da == 0.0) {
    TConf = Ta;
    return 2*Ta*A/InvLevelDensity - G4StatMFParameters::DBetaDT(Ta)*A13*A13;
  } else if (Da < 0.0) {
    do {
      Tb -= 0.5*Tb;
      ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tb);
      Db = (U+__FreeInternalE0-ECompoundNucleus)/U;
    } while (Db < 0.0);
  } else {
    do {
      Tb += 0.5*Tb;
      ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tb);
      Db = (U+__FreeInternalE0-ECompoundNucleus)/U;
    } while (Db > 0.0);
  }
  
  
  G4double eps = 1.0e-14 * abs(Tb-Ta);
  
  for (G4int i = 0; i < 1000; i++) {
    G4double Tc = (Ta+Tb)/2.0;
    if (abs(Ta-Tb) <= eps) {
      TConf = Tc;
      return 2*Tc*A/InvLevelDensity - G4StatMFParameters::DBetaDT(Tc)*A13*A13;
    }
    ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tc);
    G4double Dc = (U+__FreeInternalE0-ECompoundNucleus)/U;
    
    if (Dc == 0.0) {
      TConf = Tc;
      return 2*Tc*A/InvLevelDensity - G4StatMFParameters::DBetaDT(Tc)*A13*A13;
    }
    
    if (Da*Dc < 0.0) {
      Tb = Tc;
      Db = Dc;
    } else {
      Ta = Tc;
      Da = Dc;
    } 
  }
  
  G4cerr << 
    "G4StatMFMicrocanoncal::CalcEntropyOfCompoundNucleus: I can't calculate the temperature" 
	 << G4endl;
  
  return 0.0;
}




G4StatMFChannel *  G4StatMFMicroCanonical::ChooseAandZ(const G4Fragment & theFragment)
  // Choice of fragment atomic numbers and charges 
{
    // We choose a multiplicity (1,2,3,...) and then a channel
    G4double RandNumber = G4UniformRand();

    if (RandNumber < _WCompoundNucleus) { 
	
	G4StatMFChannel * aChannel = new G4StatMFChannel;
	aChannel->CreateFragment(theFragment.GetA(),theFragment.GetZ());
	return aChannel;
	
    } else {
	
	G4double AccumWeight = _WCompoundNucleus;
	G4std::vector<G4StatMFMicroManager*>::iterator it;
	for (it = _ThePartitionManagerVector.begin(); it != _ThePartitionManagerVector.end(); ++it) {
	    AccumWeight += (*it)->GetProbability();
	    if (RandNumber < AccumWeight) {
		return (*it)->ChooseChannel(theFragment.GetA(),theFragment.GetZ(),__MeanTemperature);
	    }
	}
	G4Exception("G4StatMFMicroCanonical::ChooseAandZ: wrong normalization!");
    }

    return 0;	
}




G4double G4StatMFMicroCanonical::CalcInvLevelDensity(const G4int anA)
{
    // Calculate Inverse Density Level
    // Epsilon0*(1 + 3 /(Af - 1))
    if (anA == 1) return 0.0;
    else return
	     G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(anA - 1.0));
}
