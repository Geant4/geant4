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
// $Id: G4HyperNucleiProperties.cc 94091 2015-11-05 15:13:52Z gcosmo $
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
// Hyper Nuclei properties based on CHIPS model (Mikhail KOSOV)
// Migrate into particles category by H.Kurashige (Sep. 2007)
// 
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HyperNucleiProperties.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"


G4double G4HyperNucleiProperties::GetNuclearMass(G4int A, G4int Z, G4int LL)
{
  if (LL==0) return  G4NucleiProperties::GetNuclearMass(A, Z);
  
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  if (A < 2 || Z < 0 || Z > A-LL || LL>A ) {
#ifdef G4VERBOSE
    if (pTable->GetVerboseLevel()>0) {
      G4cout << "G4HyperNucleiProperties::GetNuclearMass: "
	     << " Wrong values for A = " << A 
	     << " Z = " << Z 
	     << " L = " << LL
	     << G4endl;
    }
#endif
    return 0.0;
  } else if( A==2 ) {
#ifdef G4VERBOSE
    if (pTable->GetVerboseLevel()>0) {
      G4cout << "G4HyperNucleiProperties::GetNuclearMass: "
	     << " No boud state for A = " << A 
	     << " Z = " << Z 
	     << " L = " << LL
	     << G4endl;
    }
#endif
    return 0.0;
  }

  
  G4ParticleDefinition* lambda = pTable->FindParticle("lambda");
  if (lambda ==0) {
#ifdef G4VERBOSE
    if (pTable->GetVerboseLevel()>0) {
      G4cout << "G4HyperNucleiProperties::GetNuclearMass: "
	     << " Lambda is not defined " << G4endl;
    }
#endif
    return 0.0;
  }
  const G4double mLL= lambda->GetPDGMass(); // mLambda

  static const G4double b7=25.*MeV;
  static const G4double b8=10.5; // Slope
  static const G4double a2=0.13*MeV; // BindingEnergy for d+Lambda(MeV)
  static const G4double a3=2.2*MeV;  // BindingEnergy for (t/He3)+Lamb(MeV)
  static const G4double eps =0.0001*MeV; // security value (MeV)

  G4double mass =  G4NucleiProperties::GetNuclearMass(A-LL, Z); 
  // A non-"strange" nucleus
  G4double bs=0.;
  if     (A-LL ==2) bs=a2;         // for nnL,npL,ppL
  else if(A-LL ==3) bs=a3;         // for 3nL,2npL,n2pL,3pL
  else if(A-LL >3)  bs=b7*std::exp(-b8/(A-LL+1.));
  mass += LL*(mLL-bs) + eps;

  return mass;
}


G4double G4HyperNucleiProperties::GetAtomicMass(G4int A, G4int Z, G4int LL)
{
  if (A < 1 || Z < 0 || Z > A-LL || LL > A || LL <0 ) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4HyperNucleiProperties::GetAtomicMass: " 
	     << " Wrong values for A = "  << A 
	     << "  Z = " << Z 
	     << "  L = " << LL  << G4endl;
    }
#endif
    return 0.0;
    
  } else {
    G4double nuclearMass = GetNuclearMass(A, Z, LL);
    return nuclearMass + Z*electron_mass_c2
         - 1.433e-5*MeV*std::pow(G4double(Z),2.39);
  }
}

