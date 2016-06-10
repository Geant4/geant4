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
// $Id: G4StatMFMacroTriNucleon.cc 68724 2013-04-05 09:26:32Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroTriNucleon.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4StatMFMacroTriNucleon::G4StatMFMacroTriNucleon() 
  : G4VStatMFMacroCluster(3) 
{}

G4StatMFMacroTriNucleon::~G4StatMFMacroTriNucleon() 
{}

G4double 
G4StatMFMacroTriNucleon::CalcMeanMultiplicity(const G4double FreeVol, 
					      const G4double mu, 
					      const G4double nu, 
					      const G4double T)
{
  G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
  G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
  static const G4double degeneracy = 2.0+2.0;  // H3 + He3
  G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
    (1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));

  // old value was 9.224*MeV
  G4double BindingE = G4NucleiProperties::GetBindingEnergy(theA,1); 
  //	  + G4NucleiProperties::GetBindingEnergy(theA,2);

  G4double exponent = (BindingE+ theA*(mu+nu*theZARatio) - 
		       Coulomb*theZARatio*theZARatio
		       *std::pow(static_cast<G4double>(theA),5./3.))/T;
  if (exponent > 700.0) exponent = 700.0;

  _MeanMultiplicity = (degeneracy*FreeVol*theA*
		       std::sqrt(static_cast<G4double>(theA))/lambda3)*
    std::exp(exponent);
			 
  return _MeanMultiplicity;
}

G4double G4StatMFMacroTriNucleon::CalcEnergy(const G4double T)
{
  G4double Coulomb = (3./5.)*(elm_coupling/G4StatMFParameters::Getr0())*
    (1.0 - 1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.));
									
  return _Energy  = -G4NucleiProperties::GetBindingEnergy(theA,1) + 
    Coulomb * theZARatio * theZARatio 
    * std::pow(static_cast<G4double>(theA),5./3.) + (3./2.) * T;
}

G4double 
G4StatMFMacroTriNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
  G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
  G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

  G4double Entropy = 0.0;
  if (_MeanMultiplicity > 0.0)
    Entropy = _MeanMultiplicity*(2.5 + std::log((4*theA)*
        std::sqrt(static_cast<G4double>(theA))*FreeVol
        /(lambda3*_MeanMultiplicity)));
								
  return Entropy;
}
