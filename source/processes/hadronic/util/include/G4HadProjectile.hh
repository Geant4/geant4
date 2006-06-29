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
#ifndef G4HadProjectile_hh
#define G4HadProjectile_hh

#include "G4Material.hh"
class G4Track;
class G4DynamicParticle;
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
class G4ParticleDefinition;

class G4HadProjectile
{
  public:
    G4HadProjectile(const G4Track &aT);
    G4HadProjectile(const G4DynamicParticle &aT);
    const G4Material * GetMaterial() const;
    const G4ParticleDefinition * GetDefinition() const;
    const G4LorentzVector & Get4Momentum() const {return theMom;}
    G4LorentzRotation & GetTrafoToLab() {return toLabFrame;}
    G4double GetKineticEnergy() const;
    G4double GetTotalEnergy() const;
    G4double GetTotalMomentum() const;
    G4double GetGlobalTime() const {return theTime;}
    
    
  private:
  
  const G4Material * theMat;
  G4LorentzVector theOrgMom;
  G4LorentzVector theMom;
  const G4ParticleDefinition * theDef;
  G4LorentzRotation toLabFrame;
  G4double theTime;
};

#endif
