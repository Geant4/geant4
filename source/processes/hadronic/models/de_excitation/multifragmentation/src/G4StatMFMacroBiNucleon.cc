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
// $Id: G4StatMFMacroBiNucleon.cc 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroBiNucleon.hh"
#include "G4StatMFParameters.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

// Operators

static const G4double degeneracy = 3.0;

G4StatMFMacroBiNucleon & G4StatMFMacroBiNucleon::
operator=(const G4StatMFMacroBiNucleon & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroBiNucleon::operator= meant to not be accessable");
    return *this;
}

G4bool G4StatMFMacroBiNucleon::operator==(const G4StatMFMacroBiNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroBiNucleon::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFMacroBiNucleon::operator!=(const G4StatMFMacroBiNucleon & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroBiNucleon::operator!= meant to not be accessable");
    return true;
}

G4double G4StatMFMacroBiNucleon::CalcMeanMultiplicity(const G4double FreeVol, 
						      const G4double mu, 
						      const G4double nu, 
						      const G4double T)
{
  G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
  G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
    
  const G4double BindingE = G4NucleiProperties::GetBindingEnergy(theA,1); 
  //old value was 2.796*MeV
  G4double exponent = (BindingE + theA*(mu+nu*theZARatio) - 
		       G4StatMFParameters::GetCoulomb()*theZARatio*theZARatio*theA
		       *G4Pow::GetInstance()->Z23(theA))/T;

  // To avoid numerical problems
  if (exponent < -300.0) exponent = -300.0;
  else if (exponent > 300.0) exponent = 300.0;

  _MeanMultiplicity = (degeneracy*FreeVol*theA*std::sqrt((G4double)theA)/lambda3)*
    G4Exp(exponent);
			 
  return _MeanMultiplicity;
}

G4double G4StatMFMacroBiNucleon::CalcEnergy(const G4double T)
{
  _Energy  = -G4NucleiProperties::GetBindingEnergy(theA,1) + 
    G4StatMFParameters::GetCoulomb() * theZARatio * theZARatio 
    * theA*G4Pow::GetInstance()->Z23(theA) + 1.5*T;
							
  return _Energy;				
}

G4double G4StatMFMacroBiNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
  G4double Entropy = 0.0;
  if (_MeanMultiplicity > 0.0) {
    G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
    G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
    // Is this formula correct?
    Entropy = _MeanMultiplicity*(2.5+G4Log(3.0*theA*std::sqrt((G4double)theA)*FreeVol
					     /(lambda3*_MeanMultiplicity)));
  }
  return Entropy;
}
