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
