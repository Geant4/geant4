// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPContEnergyAngular.hh"

G4ReactionProduct * G4NeutronHPContEnergyAngular::Sample(G4double anEnergy, G4double massCode, G4double mass)
{
   G4ReactionProduct * result;
   G4int i,ii,iii;
   G4int it;
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
   G4double result;
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
