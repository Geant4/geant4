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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPProduct.hh" 
#include "G4Poisson.hh"
#include "G4Proton.hh"

G4ReactionProductVector * G4NeutronHPProduct::Sample(G4double anEnergy)
{
  if(theDist == NULL) return NULL;
  G4ReactionProductVector * result = new G4ReactionProductVector;
  G4double mean = theYield.GetY(anEnergy);
  G4int multi;
  multi = G4int(mean+0.0001);
  if(theMassCode==0) multi = G4Poisson(mean); // @@@@gammas. please X-check this
  theDist->SetTarget(theTarget);
  theDist->SetNeutron(theNeutron);
  G4int i;
//  G4double eMax = GetTarget()->GetMass()+GetNeutron()->GetMass()
//                  - theActualStateQValue;
  theCurrentMultiplicity = static_cast<G4int>(mean);
  G4ReactionProduct * tmp;
  for(i=0;i<multi;i++)
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    if(tmp != NULL) result->push_back(tmp);
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
    if(tmp != NULL) result->push_back(tmp);
  }
  return result;
}
