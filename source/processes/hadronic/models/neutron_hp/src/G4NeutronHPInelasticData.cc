// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPInelasticData.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4NeutronHPData.hh"

G4bool G4NeutronHPInelasticData::IsApplicable(const G4DynamicParticle*aP, const G4Element*anE)
{
  G4bool result = true;
  G4double eKin = aP->GetKineticEnergy();
  if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
  return result;
}

G4NeutronHPInelasticData::G4NeutronHPInelasticData()
{
  BuildPhysicsTable(*G4Neutron::Neutron());
}
   
G4NeutronHPInelasticData::~G4NeutronHPInelasticData()
{
  delete theCrossSections;
}
   
void G4NeutronHPInelasticData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     G4Exception("Attempt to use NeutronHP data for particles other than neutrons!!!");  
  size_t numberOfElements = G4Element::GetNumberOfElements();
  theCrossSections = new G4PhysicsTable( numberOfElements );

  // make a PhysicsVector for each element

  static const G4ElementTable *theElementTable = G4Element::GetElementTable();
  for( size_t i=0; i<numberOfElements; ++i )
  {
     G4PhysicsVector* physVec = G4NeutronHPData::
      Instance()->MakePhysicsVector((*theElementTable)[i], this);
     theCrossSections->push_back(physVec);
  }
}

void G4NeutronHPInelasticData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     G4Exception("Attempt to use NeutronHP data for particles other than neutrons!!!");  
  G4cout << "G4NeutronHPInelasticData::DumpPhysicsTable still to be implemented"<<G4endl;
}

G4double G4NeutronHPInelasticData::GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE)
{
  G4double result;
  G4bool outOfRange;
  G4int index = anE->GetIndex();
    
  result = (*((*theCrossSections)(index))).GetValue(
                             aP->GetKineticEnergy(), outOfRange);
  return result;
}
