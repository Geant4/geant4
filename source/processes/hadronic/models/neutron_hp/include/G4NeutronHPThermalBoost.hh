#ifndef G4NeutronHPThermalBoost_h
#define G4NeutronHPThermalBoost_h

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"

class G4NeutronHPThermalBoost
{
public: 
  G4double GetThermalEnergy(const G4DynamicParticle * aP, 
                            const G4Element * anE, 
			    G4double aT)
  {
    G4double theA = anE->GetN();
    G4double theZ = anE->GetZ();
    return GetThermalEnergy(aP, theA ,theZ, aT);
  }
  			    
  G4double GetThermalEnergy(const G4DynamicParticle * aP, 
                            G4double theA, G4double theZ,
			    G4double aT)
  {
    // prepare neutron
    G4double eKinetic = aP->GetKineticEnergy();
    G4ReactionProduct theNeutron( aP->GetDefinition() );
    theNeutron.SetMomentum( aP->GetMomentum() );
    theNeutron.SetKineticEnergy( eKinetic );
    G4ThreeVector neuVelo = (1./aP->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();

    // prepare properly biased thermal nucleus
    G4Nucleus aNuc;
    G4double eps = 0.0001;
    G4double eleMass; 
    eleMass = ( G4NucleiPropertiesTable::GetAtomicMass(theZ+eps, theA+eps)-
                theZ*G4Electron::ElectronDefinition()->GetPDGMass() 
  	       ) / G4Neutron::Neutron()->GetPDGMass();
  
    G4ReactionProduct aThermalNuc = aNuc.GetBiasedThermalNucleus(eleMass, neuVelo, aT);
    
    // boost to rest system and return
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, aThermalNuc);
    return boosted.GetKineticEnergy();
  }
};
#endif
