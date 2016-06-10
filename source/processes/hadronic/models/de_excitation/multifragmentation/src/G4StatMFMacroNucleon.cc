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
// $Id: G4StatMFMacroNucleon.cc 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroNucleon.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

G4StatMFMacroNucleon::G4StatMFMacroNucleon() 
  : G4VStatMFMacroCluster(1), _NeutronMeanMultiplicity(0.0),
    _ProtonMeanMultiplicity(0.0)
{}

G4StatMFMacroNucleon::~G4StatMFMacroNucleon() 
{}
	
G4double 
G4StatMFMacroNucleon::CalcMeanMultiplicity(const G4double FreeVol, 
					   const G4double mu, 
					   const G4double nu, const G4double T)
{
  if (T <= 0.0) {
    throw G4HadronicException(__FILE__, __LINE__, 
			      "G4StatMFMacroNucleon::CalcMeanMultiplicity: Temperature less or equal 0");
  }

  G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
	
  G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
	
  static const G4double degeneracy = 2.0;
	
  G4double exponent_proton = (mu + nu - G4StatMFParameters::GetCoulomb())/T;
  G4double exponent_neutron = mu/T;

  if (exponent_neutron > 300.0) exponent_neutron = 300.0;
  if (exponent_proton > 300.0) exponent_proton = 300.0;

  _NeutronMeanMultiplicity = 
    (degeneracy*FreeVol/lambda3)*G4Exp(exponent_neutron);
	
  _ProtonMeanMultiplicity = 
    (degeneracy*FreeVol/lambda3)*G4Exp(exponent_proton);

  return _MeanMultiplicity = _NeutronMeanMultiplicity + _ProtonMeanMultiplicity;
}


G4double G4StatMFMacroNucleon::CalcEnergy(const G4double T)
{
  return _Energy = G4StatMFParameters::GetCoulomb()*theZARatio*theZARatio + 1.5*T;
}

G4double 
G4StatMFMacroNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
  G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
  G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;

  G4double NeutronEntropy = 0.0;
  if (_NeutronMeanMultiplicity > 0.0)
    NeutronEntropy = _NeutronMeanMultiplicity*(2.5+G4Log(2*theA*FreeVol/
			    (lambda3*_NeutronMeanMultiplicity)));
				
  G4double ProtonEntropy = 0.0;
  if (_ProtonMeanMultiplicity > 0.0)
    ProtonEntropy = _ProtonMeanMultiplicity*(2.5+G4Log(2*theA*FreeVol/
		    (lambda3*_ProtonMeanMultiplicity)));
				
  return NeutronEntropy+ProtonEntropy;
}

