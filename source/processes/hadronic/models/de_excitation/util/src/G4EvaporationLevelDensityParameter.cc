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
// $Id: G4EvaporationLevelDensityParameter.cc 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and comment out unused part

#include "G4EvaporationLevelDensityParameter.hh"

G4EvaporationLevelDensityParameter::G4EvaporationLevelDensityParameter() 
{}

G4EvaporationLevelDensityParameter::~G4EvaporationLevelDensityParameter() 
{}

G4double 
G4EvaporationLevelDensityParameter::LevelDensityParameter(G4int A, G4int, G4double) const 
{
  return 0.1*A;
}

/*
#include "G4ShellCorrections.hh"
#include "G4SystemOfUnits.hh"

// Those values are from table 3 in 
// A.S. Iljinov et al. Nucl Phys A543 (1992) 517-557
// Table 3. alpha, beta and gamma for Cameron Shell corrections
// whithout collective effects. f-factor = 2.31.

//JMQ 17-04-08 these are not used at present in G4Evaporation 
const G4double 
G4EvaporationLevelDensityParameter::ConstEvapLevelDensityParameter = 0.125/MeV;
const G4double 
G4EvaporationLevelDensityParameter::ConstEvapLevelDensityParameter= 0.0769231/MeV;
const G4double G4EvaporationLevelDensityParameter::alpha = 0.072/MeV;
const G4double G4EvaporationLevelDensityParameter::beta = 0.257/MeV;
const G4double G4EvaporationLevelDensityParameter::gamma = 0.059/MeV;
const G4double G4EvaporationLevelDensityParameter::Bs = 1.0;

// Asymptotic Level Density Parameter
//G4double AsymptoticLDP = (alpha*A + beta*g4pow-Z23(A)*Bs)/MeV;
	
// Shape of the LDP U dependence
G4double exponent = -gamma*U;
G4double f = 1.;
if (exponent > -300.) f -= G4Exp(exponent);
G4double a = AsymptoticLDP*(1. + ShellCorrection(Z,N)*f/U);
return a;
*/
