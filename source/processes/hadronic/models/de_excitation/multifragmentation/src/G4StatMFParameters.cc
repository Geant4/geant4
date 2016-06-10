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
// $Id: G4StatMFParameters.cc 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

const G4double G4StatMFParameters::fKappa = 1.0; // dimensionless

const G4double G4StatMFParameters::fKappaCoulomb = 2.0; // dimensionless

const G4double G4StatMFParameters::fEpsilon0 = 16.0*MeV;

// Bethe-Weizsacker coefficients
const G4double G4StatMFParameters::fE0 = 16.0*MeV;

const G4double G4StatMFParameters::fBeta0 = 18.0*MeV;

const G4double G4StatMFParameters::fGamma0 = 25.0*MeV;

// Critical temperature (for liquid-gas phase transitions)
const G4double G4StatMFParameters::fCriticalTemp = 18.0*MeV;

// Nuclear radius
const G4double G4StatMFParameters::fr0 = 1.17*fermi;

const G4double G4StatMFParameters::fCoulomb = 0.6*(CLHEP::elm_coupling/fr0)*
  (1.0 - 1.0/std::pow(1.0+fKappaCoulomb,1./3.));

G4StatMFParameters::G4StatMFParameters()
{}

G4StatMFParameters::~G4StatMFParameters()
{}

G4double G4StatMFParameters::GetKappa()
{ 
  return fKappa; 
}
  
G4double G4StatMFParameters::GetKappaCoulomb()  
{ 
  return fKappaCoulomb; 
} 
  
G4double G4StatMFParameters::GetEpsilon0() 
{ 
  return fEpsilon0; 
}
  
G4double G4StatMFParameters::GetE0() 
{ 
  return fE0; 
}
  
G4double G4StatMFParameters::GetBeta0() 
{ 
  return fBeta0; 
} 
  
G4double G4StatMFParameters::GetGamma0() 
{ 
  return fGamma0; 
}
  
G4double G4StatMFParameters::GetCriticalTemp()  
{ 
  return fCriticalTemp; 
}
  
G4double G4StatMFParameters::Getr0() 
{ 
  return fr0; 
}

G4double G4StatMFParameters::GetCoulomb() 
{ 
  return fCoulomb; 
}

G4double G4StatMFParameters::Beta(G4double T) 
{
  G4double res = 0.0;
  if (T < fCriticalTemp) {
    G4double CriticalTempSqr = fCriticalTemp*fCriticalTemp;
    G4double TempSqr = T*T;
    G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
		
    res = fBeta0*tmp*std::pow(tmp,0.25);
  }
  return res;
}

G4double G4StatMFParameters::DBetaDT(G4double T) 
{
  G4double res = 0.0;
  if (T < fCriticalTemp) {
    G4double CriticalTempSqr = fCriticalTemp*fCriticalTemp;
    G4double TempSqr = T*T;
    G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
		
    res = -5.0*fBeta0*std::pow(tmp,0.25)*(CriticalTempSqr*T)/
      ((CriticalTempSqr+TempSqr)*(CriticalTempSqr+TempSqr));
  }
  return res;
}

G4double 
G4StatMFParameters::GetMaxAverageMultiplicity(G4int A) 
{
  // Maximun average multiplicity: M_0 = 2.6 for A ~ 200 
  // and M_0 = 3.3 for A <= 110
  G4double MaxAverageMultiplicity = 2.6;
  if (A <= 110) { MaxAverageMultiplicity = 3.3; }
  return MaxAverageMultiplicity;
}

