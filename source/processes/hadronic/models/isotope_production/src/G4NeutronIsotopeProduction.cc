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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4NeutronIsotopeProduction.hh"

G4NeutronIsotopeProduction::
G4NeutronIsotopeProduction()
{
  numberOfElements = G4Element::GetNumberOfElements();
  theData = new G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections> * [numberOfElements];
  for(G4int i=0; i< numberOfElements; i++)
  {
    theData[i] = new G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections>;
    if((*(G4Element::GetElementTable()))(i)->GetZ()>12 &&
       (*(G4Element::GetElementTable()))(i)->GetZ()<84) // @@@@@@ workaround to ne fixed in G4NeutronHPNames.
    {
      theData[i]->Init((*(G4Element::GetElementTable()))(i));
    }
  }
}

G4NeutronIsotopeProduction::
~G4NeutronIsotopeProduction()
{
  for(G4int i=0; i<numberOfElements; i++)
  {
    delete theData[i];
  }
  if(theData!=NULL) delete [] theData;
}


G4IsoResult * G4NeutronIsotopeProduction::
GetIsotope(const G4Track& aTrack,
           const G4Nucleus & aNucleus)
{
  // is applicable?
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4Neutron::Neutron()) return NULL;
  if(aTrack.GetKineticEnergy()>100*MeV) return NULL;

  // get the isotope
  G4Material * theMaterial = aTrack.GetMaterial();
  if(theMaterial->GetZ()<13) return NULL; //@@@@@@ workaround to ne fixed in G4NeutronHPNames.
  if(theMaterial->GetZ()>83) return NULL; //@@@@@@ workaround to ne fixed in G4NeutronHPNames.
  G4int nEleInMat = theMaterial->GetNumberOfElements();
  G4int index;
  G4double * xSec = new G4double[nEleInMat];
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
  G4IsoResult * result = theData[index]->GetProductIsotope(aTrack.GetKineticEnergy());
  return result;
}

