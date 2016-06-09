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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StatMFParameters.cc,v 1.7 2002/12/12 19:17:23 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFParameters.hh"


const G4double G4StatMFParameters::_Kappa = 1.0; // dimensionless

const G4double G4StatMFParameters::_KappaCoulomb = 2.0; // dimensionless

const G4double G4StatMFParameters::_Epsilon0 = 16.0*MeV;

// Bethe-Weizsacker coefficients
const G4double G4StatMFParameters::_E0 = 16.0*MeV;

const G4double G4StatMFParameters::_Beta0 = 18.0*MeV;

const G4double G4StatMFParameters::_Gamma0 = 25.0*MeV;

// Critical temperature (for liquid-gas phase transitions)
const G4double G4StatMFParameters::_CriticalTemp = 18.0*MeV;

// Nuclear radius
const G4double G4StatMFParameters::_r0 = 1.17*fermi;

G4double G4StatMFParameters::Beta(const G4double T)
{
  if (T > _CriticalTemp) return 0.0;
  else {
    G4double CriticalTempSqr = _CriticalTemp*_CriticalTemp;
    G4double TempSqr = T*T;
    G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
		
    return _Beta0*tmp*pow(tmp,1.0/4.0);
  }
}

G4double G4StatMFParameters::DBetaDT(const G4double T) 
{
  if (T > _CriticalTemp) return 0.0;
  else {
    G4double CriticalTempSqr = _CriticalTemp*_CriticalTemp;
    G4double TempSqr = T*T;
    G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
		
    return -5.0*_Beta0*pow(tmp,1.0/4.0)*(CriticalTempSqr*T)/
      ((CriticalTempSqr+TempSqr)*(CriticalTempSqr+TempSqr));
  }
}

G4double G4StatMFParameters::GetMaxAverageMultiplicity(const G4int A)
{
  // Maximun average multiplicity: M_0 = 2.6 for A ~ 200 
  // and M_0 = 3.3 for A <= 110
  G4double MaxAverageMultiplicity = 2.6;
  if (A <= 110) MaxAverageMultiplicity = 3.3;
  return MaxAverageMultiplicity;
}

G4StatMFParameters G4StatMFParameters::theStatMFParameters;


G4StatMFParameters * G4StatMFParameters::GetAddress()
{ return &theStatMFParameters; }

