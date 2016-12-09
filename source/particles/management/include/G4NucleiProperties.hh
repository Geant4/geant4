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
// $Id: G4NucleiProperties.hh 99159 2016-09-07 08:11:50Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
// Hadronic Process: Nuclear De-excitations by V. Lara (Oct 1998)
// Migrate into particles category by H.Kurashige (17 Nov. 98)
// Added Shell-Pairing corrections to the Cameron mass 
// excess formula by V.Lara (9 May 99)
// 
// 090331 Migrate to AME03 by Koi, Tatsumi 

#ifndef G4NucleiProperties_h
#define G4NucleiProperties_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4NucleiProperties
{
 // Class Description
 //   G4NucleiProperties is an utility class to provide mass formula of nuclei
 //   (i.e. it has static member function only)

public: 

  // Destructor
  ~G4NucleiProperties() { };

  // Default constructor ()
  G4NucleiProperties() { };


public:  // With Description

  // Give mass of nucleus A,Z
  static G4double GetNuclearMass(const G4double A, const G4double Z);
  static G4double GetNuclearMass(const G4int A, const G4int Z);

  // return 'true' if the nucleus in the stable table 
  // (i.e.in G4NucleiPropertiesTable)
  static bool IsInStableTable(const G4double A, const G4double Z);
  static bool IsInStableTable(const G4int A, const G4int Z);

  // Give binding energy 
  static G4double GetBindingEnergy(const G4int A, const G4int Z);
  static G4double GetBindingEnergy(const G4double A, const G4double Z);

  // Calculate Mass Excess of nucleus A,Z
  static G4double GetMassExcess(const G4int A, const G4int Z);
  static G4double GetMassExcess(const G4double A, const G4double Z);

  //Swich AME table in use
  static void UseOldAMETable( G4bool val = true );

private:
  // hidie methods to enforce using GetNuclearMass
  // Give mass of Atom A,Z
  static G4double GetAtomicMass(const G4double A, const G4double Z);
  
private:
  
  static G4double  AtomicMass(G4double A, G4double Z);
  
  static G4double  NuclearMass(G4double A, G4double Z);
  
  static G4double BindingEnergy(G4double A, G4double Z);
  
  static G4double MassExcess(G4double A, G4double Z);

private: 
  // table of orbit electrons mass - binding energy 
  enum  {MaxZ = 120};
  static G4ThreadLocal G4double electronMass[MaxZ];

private:
  static G4ThreadLocal G4bool   isIntialized;
  static G4ThreadLocal G4double mass_proton;
  static G4ThreadLocal G4double mass_neutron;
  static G4ThreadLocal G4double mass_deuteron;
  static G4ThreadLocal G4double mass_triton;
  static G4ThreadLocal G4double mass_alpha;
  static G4ThreadLocal G4double mass_He3;
  static G4bool use_old_evaluation;
 	
};



#endif








