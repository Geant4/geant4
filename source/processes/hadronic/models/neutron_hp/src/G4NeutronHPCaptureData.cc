// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPCaptureData.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4NeutronHPData.hh"

G4bool G4NeutronHPCaptureData::IsApplicable(const G4DynamicParticle*aP, const G4Element*anE)
{
  G4bool result = true;
  G4double eKin = aP->GetKineticEnergy();
  if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
  return result;
}

G4NeutronHPCaptureData::G4NeutronHPCaptureData()
{
  BuildPhysicsTable(*G4Neutron::Neutron());
}
   
G4NeutronHPCaptureData::~G4NeutronHPCaptureData()
{
  delete theCrossSections;
}
   
void G4NeutronHPCaptureData::BuildPhysicsTable(const G4ParticleDefinition& aP)
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

void G4NeutronHPCaptureData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     G4Exception("Attempt to use NeutronHP data for particles other than neutrons!!!");  
  G4cout << "G4NeutronHPCaptureData::DumpPhysicsTable still to be implemented"<<G4endl;
}

G4double G4NeutronHPCaptureData::GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE)
{
  G4double result;
  G4bool outOfRange;
  G4int index = anE->GetIndex();
    
  result = (*((*theCrossSections)(index))).GetValue(
                             aP->GetKineticEnergy(), outOfRange);
  return result;
}
