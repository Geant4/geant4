#include "G4NeutronIsotopeProduction.hh"

G4NeutronIsotopeProduction::
G4NeutronIsotopeProduction()
{
  numberOfElements = G4Element::GetNumberOfElements();
  theData = new G4NeutronElementIsoCrossSections * [numberOfElements];
  for(G4int i=0; i< numberOfElements; i++)
  {
    theData[i] = new G4NeutronElementIsoCrossSections;
    theData[i].Init((*(G4Element::GetElementTable()))(i));
  }
}

G4NeutronIsotopeProduction::
~G4NeutronIsotopeProduction()
{
  for(G4int i=0; i<numberOfElements; i++)
  {
    delete theData[];
  }
  if(theData!=NULL) delete theData;
}

G4double G4NeutronIsotopeProduction::
PostStepGetPhysicalInteractionLength(const G4Track& track,
			             G4double   previousStepSize,
		                     G4ForceCondition* condition)
{
  if(track.GetDynamicParticle()->GetDefinition() == G4Neutron::Neutron())
  // @@@@@@@@ take energy range of cross-sections into account.
  {
    condition = Forced;
  }
  return DBL_MAX;
}

G4VParticleChange* G4NeutronIsotopeProduction::
PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  theParticleChange.Initialize(aTrack);
  G4double nEleInMat = aTrack.GetMaterial->GetNumberOfElements();
  G4int index;
  G4double xSec = new G4double[n];
  G4double sum = 0;
  {
  for(G4int i=0; i<nEleInMat; i++)
  {
    index = theMaterial->GetElement(i)->GetIndex();
    xSec[i] = theData[index]->GetCrossSection(aTrack.GetKineticEnergy());
    sum += xSec[i];
  }
  }
  G4double random = G4UniformRand();
  G4double running = 0;
  {
  for(G4int i=0; i<nEleInMat; i++)
  {
    running += xSec[i];
    index = theMaterial->GetElement(i)->GetIndex();
    if(random<=running/sum) break;
  }
  }
  delete [] xSec;
  G4String result = theData[index]->GetProductIsotope(aTrack.GetKineticEnergy());
  theParticleChange.SetProductName(result);
  return &theParticleChange;
}

