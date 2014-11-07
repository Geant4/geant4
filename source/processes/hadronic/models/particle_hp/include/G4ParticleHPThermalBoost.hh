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
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPThermalBoost_h
#define G4ParticleHPThermalBoost_h

#include "G4HadProjectile.hh"
#include "G4Element.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"

class G4ParticleHPThermalBoost
{
public: 
  G4double GetThermalEnergy(const G4HadProjectile & aP, 
                            const G4Element * anE, 
			    G4double aT)
  {
    G4double theA = anE->GetN();
    G4double theZ = anE->GetZ();
    return GetThermalEnergy(aP, theA ,theZ, aT);
  }
  			    
  G4double GetThermalEnergy(const G4HadProjectile & aP, 
                            G4double theA, G4double theZ,
			    G4double aT)
  {
    // prepare neutron
    G4double eKinetic = aP.GetKineticEnergy();
    G4ReactionProduct theNeutronRP( const_cast<G4ParticleDefinition *>(aP.GetDefinition()) );
    theNeutronRP.SetMomentum( aP.Get4Momentum().vect() );
    theNeutronRP.SetKineticEnergy( eKinetic );
    G4ThreeVector neuVelo = (1./aP.GetDefinition()->GetPDGMass())*theNeutronRP.GetMomentum();

    // prepare properly biased thermal nucleus
    G4Nucleus aNuc;
    G4double eps = 0.0001;
    G4double eleMass; 
    eleMass = ( G4NucleiProperties::GetNuclearMass( static_cast<G4int>(theA+eps) , static_cast<G4int>(theZ+eps) ) ) / G4Neutron::Neutron()->GetPDGMass();
  
    G4ReactionProduct aThermalNuc = aNuc.GetBiasedThermalNucleus(eleMass, neuVelo, aT);
    
    // boost to rest system and return
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutronRP, aThermalNuc);
    return boosted.GetKineticEnergy();
  }
};
#endif
