// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPFissionData.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4NeutronHPData.hh"

G4bool G4NeutronHPFissionData::IsApplicable(const G4DynamicParticle*aP, const G4Element*anE)
{
  G4bool result = true;
  G4double eKin = aP->GetKineticEnergy();
  if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
  return result;
}

G4NeutronHPFissionData::G4NeutronHPFissionData()
{
  theCrossSections = NULL;
  BuildPhysicsTable(*G4Neutron::Neutron());
}
   
G4NeutronHPFissionData::~G4NeutronHPFissionData()
{
  if(theCrossSections!=NULL) delete theCrossSections;
}
   
void G4NeutronHPFissionData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     G4Exception("Attempt to use NeutronHP data for particles other than neutrons!!!");  
  G4int numberOfElements = G4Element::GetNumberOfElements();
  theCrossSections = new G4PhysicsTable( numberOfElements );

  // make a PhysicsVector for each element

  static const G4ElementTable *theElementTable = G4Element::GetElementTable();
  for( G4int i=0; i<numberOfElements; ++i )
    (*theCrossSections)(i) =
      G4NeutronHPData::
      Instance()->MakePhysicsVector((*theElementTable)[i], this);
}

void G4NeutronHPFissionData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     G4Exception("Attempt to use NeutronHP data for particles other than neutrons!!!");  
  G4cout << "G4NeutronHPFissionData::DumpPhysicsTable still to be implemented"<<G4endl;
}

G4double G4NeutronHPFissionData::GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE)
{
  G4double result;
  G4bool outOfRange;
  G4int index = anE->GetIndex();
    
  if(anE->GetZ()<90) return 0;
  result = (*((*theCrossSections)(index))).GetValue(
                             aP->GetTotalMomentum()/GeV, outOfRange);
  return result;
}
