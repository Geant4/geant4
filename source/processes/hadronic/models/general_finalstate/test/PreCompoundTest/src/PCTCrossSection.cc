#include "PCTCrossSection.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4DynamicParticle.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

#include "PCTTarget.hh"
#include "PCTProjectile.hh"

G4double PCTCrossSection::CalculateCrossSection(PCTProjectile * projectile, PCTTarget * target)
{
    if (target->IsValid())
    {
	if (target->IsNucleus()) this->CalculateNucleusXS(projectile,target);
	else if (target->IsNatural()) this->CalculateElementXS(projectile,target);
	else CrossSection = -1.0;
    }
    return this->GetCrossSection();
}

void PCTCrossSection::CalculateNucleusXS(PCTProjectile * projectile, PCTTarget * target)
{
    // Calculate the Cross Section
    G4int A0 = target->GetA();
    G4int Z0 = target->GetZ();

    G4double massofmole = G4ParticleTable::GetParticleTable()->GetIonTable()->
	GetIonMass(Z0,A0) * ((Avogadro*mole)/(g*c_squared))*(g/mole); // g/mole  

    
    G4Isotope * isotope = new G4Isotope("Isotope", Z0, A0, massofmole );
    G4int ncomponents = 1;
    G4Element * element = new G4Element("Elelement", "El", ncomponents);
    element->AddIsotope(isotope, 100*perCent);
    G4DynamicParticle * aDP = new G4DynamicParticle(projectile->GetDefinition(), 
						    projectile->GetMomentum());
    G4VCrossSectionDataSet * aXS = this->DataSetAdapter(projectile->GetDefinition());
    aXS->BuildPhysicsTable(*(projectile->GetDefinition()));
    CrossSection = aXS->GetCrossSection(aDP,element);
    delete aXS;
    delete element;
    delete isotope;
    return;
}

void PCTCrossSection::CalculateElementXS(PCTProjectile * projectile, PCTTarget * target)
{
    // Calculate the Cross Section
    G4int A0 = target->GetA();
    G4int Z0 = target->GetZ();

    G4double massofmole = G4ParticleTable::GetParticleTable()->GetIonTable()->
	GetIonMass(Z0,A0) * ((Avogadro*mole)/(g*c_squared))*(g/mole); // g/mole  

    
    G4Isotope * isotope = new G4Isotope("Isotope", Z0, A0, massofmole );
    G4int ncomponents = 1;
    G4Element * element = new G4Element("Elelement", "El", ncomponents);
    element->AddIsotope(isotope, 100*perCent);
    G4DynamicParticle * aDP = new G4DynamicParticle(projectile->GetDefinition(), 
						    projectile->GetMomentum());
    G4VCrossSectionDataSet * aXS = this->DataSetAdapter(projectile->GetDefinition());
    aXS->BuildPhysicsTable(*(projectile->GetDefinition()));
    CrossSection = aXS->GetCrossSection(aDP,element);
    delete aXS;
    delete element;
    delete isotope;
    return;
}


G4VCrossSectionDataSet * PCTCrossSection::DataSetAdapter(const G4ParticleDefinition * particle)
{
    G4String name = particle->GetParticleName();
    if (name == "proton")
	return new G4ProtonInelasticCrossSection();
    else if (name == "neutron")
	return new G4NeutronInelasticCrossSection();
    else
	return 0;
}
