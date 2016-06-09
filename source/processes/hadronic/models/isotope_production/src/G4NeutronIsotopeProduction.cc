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
#include "G4NeutronIsotopeProduction.hh"

G4NeutronIsotopeProduction::
G4NeutronIsotopeProduction()
{
  numberOfElements = G4Element::GetNumberOfElements();
  theData = new G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections> * [numberOfElements];
  for(G4int i=0; i< numberOfElements; i++)
  {
    theData[i] = new G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections>;
    if((*(G4Element::GetElementTable()))[i]->GetZ()>12 &&
       (*(G4Element::GetElementTable()))[i]->GetZ()<84) // @@@@@@ workaround to ne fixed in G4NeutronHPNames.
    {
      theData[i]->Init((*(G4Element::GetElementTable()))[i]);
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
           const G4Nucleus & )
{
  // is applicable?
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4Neutron::Neutron()) return NULL;
  if(aTrack.GetKineticEnergy()>100*MeV) return NULL;

  // get the isotope
  G4Material * theMaterial = aTrack.GetMaterial();
  G4int nEleInMat = theMaterial->GetNumberOfElements();
  for(G4int check = 0; check<nEleInMat; check++)
  {
    if(theMaterial->GetElement(check)->GetZ()<13) return NULL; //@@@@@@ workaround to ne fixed in G4NeutronHPNames.
    if(theMaterial->GetElement(check)->GetZ()>83) return NULL; //@@@@@@ workaround to ne fixed in G4NeutronHPNames.
  }
  G4int index(0);
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

