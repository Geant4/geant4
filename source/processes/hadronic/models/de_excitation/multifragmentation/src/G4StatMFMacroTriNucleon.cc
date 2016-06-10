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
// $Id: G4StatMFMacroTriNucleon.cc 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroTriNucleon.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

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
  static const G4double degeneracy = 4.0;  // H3 + He3

  // old value was 9.224*MeV
  G4double BindingE = G4NucleiProperties::GetBindingEnergy(theA,1); 
  //	  + G4NucleiProperties::GetBindingEnergy(theA,2);

  G4double exponent = (BindingE+ theA*(mu+nu*theZARatio) - 
		       G4StatMFParameters::GetCoulomb()*theZARatio*theZARatio
		       *theA*G4Pow::GetInstance()->Z23(theA))/T;
  if (exponent > 300.0) exponent = 300.0;

  _MeanMultiplicity = (degeneracy*FreeVol*theA*std::sqrt((G4double)theA)/lambda3)*
    G4Exp(exponent);
			 
  return _MeanMultiplicity;
}

G4double G4StatMFMacroTriNucleon::CalcEnergy(const G4double T)
{
  return _Energy  = -G4NucleiProperties::GetBindingEnergy(theA,1) + 
    G4StatMFParameters::GetCoulomb() * theZARatio * theZARatio 
    * theA*G4Pow::GetInstance()->Z23(theA) + 1.5 * T;
}

G4double 
G4StatMFMacroTriNucleon::CalcEntropy(const G4double T, const G4double FreeVol)
{
  G4double Entropy = 0.0;
  if (_MeanMultiplicity > 0.0) {
    G4double ThermalWaveLenght = 16.15*fermi/std::sqrt(T);
    G4double lambda3 = ThermalWaveLenght*ThermalWaveLenght*ThermalWaveLenght;
    Entropy = _MeanMultiplicity*(2.5 + G4Log((4*theA)*std::sqrt((G4double)theA)
					     *FreeVol/(lambda3*_MeanMultiplicity)));
  }				
  return Entropy;
}
