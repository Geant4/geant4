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
// $Id: G4NucleiProperties.cc 99159 2016-09-07 08:11:50Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
// Migrate into particles category by H.Kurashige (17 Nov. 98)
// Added Shell-Pairing corrections to the Cameron mass 
// excess formula by V.Lara (9 May 99)
// 090331 Migrate to AME03 by Koi, Tatsumi 

#include "G4NucleiProperties.hh"

#include "G4NucleiPropertiesTableAME03.hh"
#include "G4NucleiPropertiesTableAME12.hh"
#include "G4NucleiPropertiesTheoreticalTable.hh"
#include "G4ParticleTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4double G4NucleiProperties::mass_proton = -1.;
G4ThreadLocal G4double G4NucleiProperties::mass_neutron = -1.;
G4ThreadLocal G4double G4NucleiProperties::mass_deuteron = -1.;
G4ThreadLocal G4double G4NucleiProperties::mass_triton = -1.;
G4ThreadLocal G4double G4NucleiProperties::mass_alpha = -1.;
G4ThreadLocal G4double G4NucleiProperties::mass_He3 = -1.;
#ifndef G4NucleiProperties_USE_OLD_AME_TABLE
G4bool G4NucleiProperties::use_old_evaluation = false;
#else
G4bool G4NucleiProperties::use_old_evaluation = true;
#endif

G4double G4NucleiProperties::GetNuclearMass(const G4double A, const G4double Z)
{
  G4double mass=0.0;

  if (std::fabs(A - G4int(A)) > 1.e-10) {
    mass = NuclearMass(A,Z);
 
  } else {
    // use mass table
    G4int iZ = G4int(Z);
    G4int iA = G4int(A);
    mass =GetNuclearMass(iA,iZ);
  }
  
   return mass;
}


G4double G4NucleiProperties::GetNuclearMass(const G4int A, const G4int Z)
{
  if (mass_proton  <= 0.0 ) {
    const G4ParticleDefinition * nucleus = 0;
    nucleus = G4ParticleTable::GetParticleTable()->FindParticle("proton"); // proton 
    if (nucleus!=0) mass_proton = nucleus->GetPDGMass();
    nucleus = G4ParticleTable::GetParticleTable()->FindParticle("neutron"); // neutron 
    if (nucleus!=0) mass_neutron = nucleus->GetPDGMass();
    nucleus = G4ParticleTable::GetParticleTable()->FindParticle("deuteron"); // deuteron 
    if (nucleus!=0) mass_deuteron = nucleus->GetPDGMass();
    nucleus = G4ParticleTable::GetParticleTable()->FindParticle("triton"); // triton 
    if (nucleus!=0) mass_triton = nucleus->GetPDGMass();
    nucleus = G4ParticleTable::GetParticleTable()->FindParticle("alpha"); // alpha 
    if (nucleus!=0) mass_alpha = nucleus->GetPDGMass();
    nucleus = G4ParticleTable::GetParticleTable()->FindParticle("He3"); // He3 
    if (nucleus!=0) mass_He3 = nucleus->GetPDGMass();

  }

  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4NucleiProperties::GetNuclearMass: Wrong values for A = " << A 
	     << " and Z = " << Z << G4endl;
    }
#endif    
    return 0.0;
  }
  
  G4double mass= -1.;
  if ( (Z<=2) ) {
    // light nuclei
    if ( (Z==1)&&(A==1) ) {
      mass = mass_proton;
    } else if ( (Z==0)&&(A==1) ) {
      mass = mass_neutron;
    } else if ( (Z==1)&&(A==2) ) {
      mass = mass_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      mass = mass_triton;
    } else if ( (Z==2)&&(A==4) ) {
      mass = mass_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      mass = mass_He3;
    }
  }
  
  if (mass < 0.) {
    G4bool inAMETable = false;
    if ( ! use_old_evaluation ) {
       inAMETable = G4NucleiPropertiesTableAME12::IsInTable(Z,A);
    } else {
       inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    }
    if ( inAMETable ) {
      // AME table
      if ( ! use_old_evaluation ) {
         mass = G4NucleiPropertiesTableAME12::GetNuclearMass(Z,A);
      } else {
         mass = G4NucleiPropertiesTableAME03::GetNuclearMass(Z,A);
      }
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A)){
      // Theoretical table
      mass = G4NucleiPropertiesTheoreticalTable::GetNuclearMass(Z,A);
    } else {
      mass = NuclearMass(G4double(A),G4double(Z));
    }
  }

  if (mass < 0.) mass = 0.0;
  return mass;
}

G4bool G4NucleiProperties::IsInStableTable(const G4double A, const G4double Z)
{
  G4int iA = G4int(A);
  G4int iZ = G4int(Z);
  return IsInStableTable(iA, iZ);
}

G4bool G4NucleiProperties::IsInStableTable(const G4int A, const int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4NucleiProperties::IsInStableTable: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;	
    }
#endif 
    return false;
  } 

  if ( ! use_old_evaluation ) {
    return G4NucleiPropertiesTableAME12::IsInTable(Z,A);
  } else {
    return G4NucleiPropertiesTableAME03::IsInTable(Z,A);
  }
}

G4double G4NucleiProperties::GetMassExcess(const G4double A, const G4double Z)
{
  G4int iA = G4int(A);
  G4int iZ = G4int(Z);
  return GetMassExcess(iA,iZ);
}

G4double G4NucleiProperties::GetMassExcess(const G4int A, const G4int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4NucleiProperties::GetMassExccess: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;
    }
#endif    
    return 0.0;
    
  } else {

    G4bool inAMETable = false;
    if ( ! use_old_evaluation ) {
       inAMETable = G4NucleiPropertiesTableAME12::IsInTable(Z,A);
    } else {
       inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    }
    if (inAMETable){
      // AME table
      if ( ! use_old_evaluation ) {
         return G4NucleiPropertiesTableAME12::GetMassExcess(Z,A);
      } else {
         return G4NucleiPropertiesTableAME03::GetMassExcess(Z,A);
      }
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A)){
      return G4NucleiPropertiesTheoreticalTable::GetMassExcess(Z,A);
    } else {
      return MassExcess(A,Z);
    }
  }

}


G4double G4NucleiProperties::GetAtomicMass(const G4double A, const G4double Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4NucleiProperties::GetAtomicMass: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;	
    }
#endif 
    return 0.0;

  } else if (std::fabs(A - G4int(A)) > 1.e-10) {
    return AtomicMass(A,Z);

  } else {
    G4int iA = G4int(A);
    G4int iZ = G4int(Z);
    G4bool inAMETable = false;
    if ( ! use_old_evaluation ) {
       inAMETable = G4NucleiPropertiesTableAME12::IsInTable(Z,A);
    } else {
       inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    }
    if (inAMETable) {
      if ( ! use_old_evaluation ) {
        return G4NucleiPropertiesTableAME12::GetAtomicMass(Z,A);
      } else {
        return G4NucleiPropertiesTableAME03::GetAtomicMass(Z,A);
      }
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(iZ,iA)){
      return G4NucleiPropertiesTheoreticalTable::GetAtomicMass(iZ,iA);
    } else {
      return AtomicMass(A,Z);
    }
  }
}

G4double G4NucleiProperties::GetBindingEnergy(const G4double A, const G4double Z)
{
  G4int iA = G4int(A);
  G4int iZ = G4int(Z);
  return GetBindingEnergy(iA,iZ);
}

G4double G4NucleiProperties::GetBindingEnergy(const G4int A, const G4int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4NucleiProperties::GetMassExccess: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;
    }
#endif
    return 0.0;

  } else {
    G4bool inAMETable = false;
    if ( ! use_old_evaluation ) {
       inAMETable = G4NucleiPropertiesTableAME12::IsInTable(Z,A);
    } else {
       inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    }
    if (inAMETable) {
      if ( ! use_old_evaluation ) {
        return G4NucleiPropertiesTableAME12::GetBindingEnergy(Z,A);
      } else {
        return G4NucleiPropertiesTableAME03::GetBindingEnergy(Z,A);
      }
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A)) {
      return G4NucleiPropertiesTheoreticalTable::GetBindingEnergy(Z,A);
    }else {
      return BindingEnergy(A,Z);
    }

  }
}


G4double G4NucleiProperties::MassExcess(G4double A, G4double Z) 
{
  return GetAtomicMass(A,Z) - A*amu_c2;
}

G4double  G4NucleiProperties::AtomicMass(G4double A, G4double Z)
{
  //const G4double hydrogen_mass_excess;
  //const G4double neutron_mass_excess;  
  G4double hydrogen_mass_excess;
  G4double neutron_mass_excess;  
  if ( ! use_old_evaluation ) {
    hydrogen_mass_excess = G4NucleiPropertiesTableAME12::GetMassExcess(1,1);
    neutron_mass_excess = G4NucleiPropertiesTableAME12::GetMassExcess(0,1);
  } else {
    hydrogen_mass_excess = G4NucleiPropertiesTableAME03::GetMassExcess(1,1);
    neutron_mass_excess = G4NucleiPropertiesTableAME03::GetMassExcess(0,1);
  }

  G4double mass =
      (A-Z)*neutron_mass_excess + Z*hydrogen_mass_excess - BindingEnergy(A,Z) + A*amu_c2;

  return mass;
}

G4double  G4NucleiProperties::NuclearMass(G4double A, G4double Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "G4NucleiProperties::NuclearMass: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;
    }
#endif 
    return 0.0;
  }

  G4double mass = AtomicMass(A,Z);
  // atomic mass is converted to nuclear mass according formula in  AME03 and 12
  mass -= Z*electron_mass_c2;
  mass += ( 14.4381*std::pow ( Z , 2.39 ) + 1.55468*1e-6*std::pow ( Z , 5.35 ) )*eV;      

  return mass;
}

G4double  G4NucleiProperties::BindingEnergy(G4double A, G4double Z)
{ 
  //
  // Weitzsaecker's Mass formula
  //
  G4int Npairing = G4int(A-Z)%2;                  // pairing
  G4int Zpairing = G4int(Z)%2;
  G4double binding =
      - 15.67*A                           // nuclear volume
      + 17.23*std::pow(A,2./3.)                // surface energy
      + 93.15*((A/2.-Z)*(A/2.-Z))/A       // asymmetry
      + 0.6984523*Z*Z*std::pow(A,-1./3.);      // coulomb
  if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / std::sqrt(A);  // pairing

  return -binding*MeV;
}

void G4NucleiProperties::UseOldAMETable( G4bool val )
{
   use_old_evaluation = val;
}
