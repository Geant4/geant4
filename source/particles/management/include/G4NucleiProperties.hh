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
// G4NucleiProperties
//
// Class description:
//
// Utility class to provide mass formula of nuclei (i.e. it has static
// member function only).

// Author: V.Lara, October 1998
// History:
// - 17.11.1998, H.Kurashige - Migrated into particles category
// - 31.03.2009, T.Koi - Migrated to AME03
// --------------------------------------------------------------------
#ifndef G4NucleiProperties_hh
#define G4NucleiProperties_hh 1

#include "globals.hh"
#include "G4ios.hh"

class G4NucleiProperties
{
  public: 

    G4NucleiProperties() {}
   ~G4NucleiProperties() {}

    static G4double GetNuclearMass(const G4double A, const G4double Z);
    static G4double GetNuclearMass(const G4int A, const G4int Z);
      // Give mass of nucleus A,Z

    static G4bool IsInStableTable(const G4double A, const G4double Z);
    static G4bool IsInStableTable(const G4int A, const G4int Z);
      // Return 'true' if the nucleus is in the stable table 
      // (i.e. in G4NucleiPropertiesTable)

    static G4double GetBindingEnergy(const G4int A, const G4int Z);
    static G4double GetBindingEnergy(const G4double A, const G4double Z);
      // Give binding energy 

    static G4double GetMassExcess(const G4int A, const G4int Z);
    static G4double GetMassExcess(const G4double A, const G4double Z);
      // Calculate Mass Excess of nucleus A,Z

  private:

    // Hidden methods to enforce using GetNuclearMass

    static G4double GetAtomicMass(const G4double A, const G4double Z);
      // Give mass of Atom A,Z
  
    static G4double AtomicMass(G4double A, G4double Z);
  
    static G4double NuclearMass(G4double A, G4double Z);
  
    static G4double BindingEnergy(G4double A, G4double Z);
  
    static G4double MassExcess(G4double A, G4double Z);

  private: 

    enum  {MaxZ = 120};
    static G4ThreadLocal G4double electronMass[MaxZ];
      // table of orbit electrons mass - binding energy 

    static G4ThreadLocal G4bool   isIntialized;
    static G4ThreadLocal G4double mass_proton;
    static G4ThreadLocal G4double mass_neutron;
    static G4ThreadLocal G4double mass_deuteron;
    static G4ThreadLocal G4double mass_triton;
    static G4ThreadLocal G4double mass_alpha;
    static G4ThreadLocal G4double mass_He3;
};

#endif
