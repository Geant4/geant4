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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//
#ifndef G4ParticleHPNames_h
#define G4ParticleHPNames_h 1

#include "globals.hh"
#include "G4ParticleHPDataUsed.hh"

class G4ParticleHPManager;

class G4ParticleHPNames
{
  public:
    explicit G4ParticleHPNames(G4int maxOffSet = 5);
    ~G4ParticleHPNames() = default;

    G4ParticleHPDataUsed GetName(G4int A, G4int Z, const G4String& base,
                                 const G4String& rest, G4bool& active)
    {
      return GetName(A, Z, 0, base, rest, active);
    }
    G4ParticleHPDataUsed GetName(G4int A, G4int Z, G4int M,
                                 const G4String& base, const G4String& rest, G4bool& active);
    G4String GetName(G4int i);

    G4String itoa(G4int current);

    void SetMaxOffSet(G4int anOffset) { theMaxOffSet = anOffset; }

    G4ParticleHPNames(G4ParticleHPNames&) = delete;
    G4ParticleHPNames& operator=(const G4ParticleHPNames &right) = delete;

  private:
    G4ParticleHPManager* fManager;
    static const G4String theString[100];
    G4int theMaxOffSet;
};

#endif
