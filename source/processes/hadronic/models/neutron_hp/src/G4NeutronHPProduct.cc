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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080718 As for secondary photons, if its mean value has a value of integer, 
//        then a sampling of multiplicity that based on Poisson Distribution 
//        is not carried out and the mean is used as a multiplicity.
//        modified by T. Koi.
// 080721 Using ClearHistories() methodl for limiting the sum of secondary energies
//        modified by T. Koi.
// 080901 bug fix of too many secnodaries production in nd reactinos by T. Koi
//
#include "G4NeutronHPProduct.hh" 
#include "G4Poisson.hh"
#include "G4Proton.hh"

G4ReactionProductVector * G4NeutronHPProduct::Sample(G4double anEnergy)
{
  if(theDist == 0) { return 0; }
  G4ReactionProductVector * result = new G4ReactionProductVector;
  G4double mean = theYield.GetY(anEnergy);
  G4int multi;
  multi = G4int(mean+0.0001);
  //if(theMassCode==0) multi = G4Poisson(mean); // @@@@gammas. please X-check this
  //080718
  if ( theMassCode == 0 )
  { 
     if ( G4int ( mean ) == mean )
     {
        multi = (G4int) mean;
     }
     else
     {
        multi = G4Poisson ( mean ); 
     }
  }
  theDist->SetTarget(theTarget);
  theDist->SetNeutron(theNeutron);
  G4int i;
//  G4double eMax = GetTarget()->GetMass()+GetNeutron()->GetMass()
//                  - theActualStateQValue;
  theCurrentMultiplicity = static_cast<G4int>(mean);
  G4ReactionProduct * tmp;
  theDist->ClearHistories();
  for(i=0;i<multi;i++)
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    if(tmp != 0) { result->push_back(tmp); }
  }
  if(multi == 0) 
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    delete  tmp;
  }
/*
//080901 TK Comment out, too many secondaries are produced in deuteron reactions
  if(theTarget->GetMass()<2*GeV) // @@@ take care of residuals in all cases
  {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    tmp->SetDefinition(G4Proton::Proton());
    if(tmp != 0) { result->push_back(tmp); }
  }
*/
  return result;
}
