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
// $Id: G4StatMFMacroMultiNucleon.cc 100379 2016-10-19 15:05:35Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 25.07.08 I.Pshenichnov (in collaboration with Alexander Botvina and Igor 
//          Mishustin (FIAS, Frankfurt, INR, Moscow and Kurchatov Institute, 
//          Moscow, pshenich@fias.uni-frankfurt.de) fixed computation of the 
//          symmetry energy 

#include "G4StatMFMacroMultiNucleon.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

// Default constructor
G4StatMFMacroMultiNucleon::
G4StatMFMacroMultiNucleon() :
    G4VStatMFMacroCluster(0)  // Beacuse the def. constr. of base class is private
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::default_constructor meant to not be accessable");
}

// Copy constructor
G4StatMFMacroMultiNucleon::
G4StatMFMacroMultiNucleon(const G4StatMFMacroMultiNucleon & ) :
    G4VStatMFMacroCluster(0)  // Beacuse the def. constr. of base class is private
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMacroMultiNucleon & G4StatMFMacroMultiNucleon::
operator=(const G4StatMFMacroMultiNucleon & ) 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::operator= meant to not be accessable");
    return *this;
}

G4bool G4StatMFMacroMultiNucleon::operator==(const G4StatMFMacroMultiNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::operator== meant to not be accessable");
    return false;
}
 
G4bool G4StatMFMacroMultiNucleon::operator!=(const G4StatMFMacroMultiNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroMultiNucleon::operator!= meant to not be accessable");
    return true;
}

G4double G4StatMFMacroMultiNucleon::CalcMeanMultiplicity(const G4double FreeVol, 
							 const G4double mu,
							 const G4double nu, 
							 const G4double T)
{
  G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);	
  G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double A23 = g4calc->Z23(theA);
	
  G4double exponent = (mu + nu*theZARatio+ G4StatMFParameters::GetE0() 
		       + T*T/_InvLevelDensity 
		       - G4StatMFParameters::GetGamma0()*(1.0 - 2.0*theZARatio)*
		       (1.0 - 2.0*theZARatio))*theA
    - G4StatMFParameters::Beta(T)*A23 
    - G4StatMFParameters::GetCoulomb()*theZARatio*theZARatio*A23*theA;
	
  exponent /= T;
	
  if (exponent > 30.0) exponent = 30.0;
	
  _MeanMultiplicity = std::max((FreeVol * theA * std::sqrt((G4double)theA)/lambda3) *
			       G4Exp(exponent),1.0e-30);
  return _MeanMultiplicity;	
}

G4double G4StatMFMacroMultiNucleon::CalcZARatio(const G4double nu)
{
  G4double den = 8*G4StatMFParameters::GetGamma0() 
    + 2*G4StatMFParameters::GetCoulomb()*G4Pow::GetInstance()->Z23(theA);
  theZARatio = (4.0*G4StatMFParameters::GetGamma0()+nu)/den;
  return theZARatio;
}

G4double G4StatMFMacroMultiNucleon::CalcEnergy(const G4double T)
{
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double A23 = g4calc->Z23(theA);

  // Volume term 
  G4double EVol = theA * (T*T/_InvLevelDensity - G4StatMFParameters::GetE0());
	
  // Symmetry term
  G4double ESym = theA * G4StatMFParameters::GetGamma0() 
    *(1. - 2.* theZARatio) * (1. - 2.* theZARatio);
	
  // Surface term
  G4double ESurf = A23*(G4StatMFParameters::Beta(T) - T*G4StatMFParameters::DBetaDT(T));
 
  // Coulomb term
  G4double ECoul = G4StatMFParameters::GetCoulomb()*A23*theA*theZARatio*theZARatio;
	
  // Translational term
  G4double ETrans = 1.5*T;
  return _Energy = EVol + ESurf + ECoul + ETrans + ESym;
}

G4double G4StatMFMacroMultiNucleon::CalcEntropy(const G4double T, 
						const G4double FreeVol)
{
  G4double Entropy = 0.0;
  if (_MeanMultiplicity > 0.0) {

    G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
    G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
    // Volume term
    G4double SV = 2.0*theA*T/_InvLevelDensity;
		
    // Surface term
    G4double SS = -G4StatMFParameters::DBetaDT(T)*G4Pow::GetInstance()->Z23(theA);
		
    // Translational term
    G4double ST = 2.5 + G4Log(FreeVol * std::sqrt((G4double)theA) * theA
			      /(lambda3*_MeanMultiplicity));
				
    Entropy = _MeanMultiplicity*(SV + SS + ST);
  }								
  return Entropy;
}
