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
//
#include "G4NeutronHPContEnergyAngular.hh"

G4ReactionProduct * G4NeutronHPContEnergyAngular::Sample(G4double anEnergy, G4double massCode, G4double mass)
{
   G4ReactionProduct * result;
   G4int i(0);
   G4int it(0);
   for(i=0;i<nEnergy;i++)
   {
     it = i;
     if(theAngular[i].GetEnergy()>anEnergy) break;
   }
   G4double targetMass = GetTarget()->GetMass();
   if(it==0)
   {
     theAngular[0].SetTarget(GetTarget());
     theAngular[0].SetTargetCode(theTargetCode);
     theAngular[0].SetPrimary(GetNeutron());
     result = theAngular[0].Sample(anEnergy, massCode, targetMass, 
                                  theAngularRep, theInterpolation);
     currentMeanEnergy = theAngular[0].MeanEnergyOfThisInteraction();
   }
   else
   {
     // interpolation through alternating sampling. This needs improvement @@@
     // This is the cause of the He3 problem !!!!!!!!
     // See to it, if you can improve this.
     G4double random = G4UniformRand();
     G4double deltaE = theAngular[it].GetEnergy()-theAngular[it-1].GetEnergy();
     G4double offset = theAngular[it].GetEnergy()-anEnergy;
     if(random<offset/deltaE) it--; 
     theAngular[it].SetTarget(GetTarget());
     theAngular[it].SetTargetCode(theTargetCode);
     theAngular[it].SetPrimary(GetNeutron());
     result = theAngular[it].Sample(anEnergy, massCode, targetMass, 
                                    theAngularRep, theInterpolation);
     currentMeanEnergy = theAngular[it].MeanEnergyOfThisInteraction();
   }
   return result;
}

G4double G4NeutronHPContEnergyAngular::
MeanEnergyOfThisInteraction()
{
   G4double result(0);
   if(currentMeanEnergy<-1)
   {
     G4Exception("G4NeutronHPContEnergyAngular: Logical error in Product class");
   }
   else
   {
     result = currentMeanEnergy;
   }
   currentMeanEnergy = -2;
   return result;
}
