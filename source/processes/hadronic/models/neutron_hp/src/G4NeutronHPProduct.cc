// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPProduct.hh" 
#include "Randomize.hh"
#include "G4Proton.hh"

  G4NeutronHPProduct::G4NeutronHPProduct()
  {
    theDist = NULL;
  }
  
  G4NeutronHPProduct::~G4NeutronHPProduct()
  {
    if(theDist != NULL) delete theDist;
  }

G4ReactionProductVector * G4NeutronHPProduct::Sample(G4double anEnergy)
{
  if(theDist == NULL) return NULL;
  G4ReactionProductVector * result = new G4ReactionProductVector;
  G4double mean = theYield.GetY(anEnergy);
  G4int multi;
  multi = G4int(mean+0.0001);
  if(theMassCode==0) multi = RandPoisson::shoot(mean); // @@@@gammas. please X-check this
  theDist->SetTarget(theTarget);
  theDist->SetNeutron(theNeutron);
  G4int i;
  G4double eMax = GetTarget()->GetMass()+GetNeutron()->GetMass()
                  - theActualStateQValue;
  theCurrentMultiplicity = mean;
  G4ReactionProduct * tmp;
  for(i=0;i<multi;i++)
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    if(tmp != NULL) result->insert(tmp);
  }
  if(multi == 0) 
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    delete  tmp;
  }
  if(theTarget->GetMass()<2*GeV) // @@@ take care of residuals in all cases
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    tmp->SetDefinition(G4Proton::Proton());
    if(tmp != NULL) result->insert(tmp);
  }
  return result;
}
