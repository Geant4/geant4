// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPIsotropic.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"

void G4NeutronHPIsotropic::Init(ifstream & aDataFile)
{
}

G4ReactionProduct * G4NeutronHPIsotropic::Sample(G4double anEnergy, G4double massCode, G4double mass)
{
   G4ReactionProduct * result = new G4ReactionProduct;
   G4int Z = massCode/1000;
   G4int A = massCode-1000*Z;

   if(massCode==0)
   {
     result->SetDefinition(G4Gamma::Gamma());
   }
   else if(A==0)
   {
     result->SetDefinition(G4Electron::Electron());     
     if(Z==1) result->SetDefinition(G4Positron::Positron());
   }
   else if(A==1)
   {
     result->SetDefinition(G4Neutron::Neutron());
     if(Z==1) result->SetDefinition(G4Proton::Proton());
   }
   else if(A==2)
   {
     result->SetDefinition(G4Deuteron::Deuteron());      
   }
   else if(A==3)
   {
     result->SetDefinition(G4Triton::Triton());  
     if(Z==2) result->SetDefinition(G4He3::He3());
   }
   else if(A==4)
   {
     result->SetDefinition(G4Alpha::Alpha());
     if(Z!=2) G4Exception("Unknown ion case 1");    
   }
   else
   {
     G4Exception("G4NeutronHPIsotropic: Unknown ion case 2");
   }

   G4double cosTh = G4UniformRand();
   G4double phi = twopi*G4UniformRand();
   G4double theta = acos(cosTh);
   G4double sinth = sin(theta);
   
// we need the the Q value of the reaction
   result->SetKineticEnergy(max(0.001*MeV, anEnergy+GetQValue()));
   G4double mtot = result->GetTotalMomentum(); 
   G4ThreeVector tempVector(mtot*sinth*cos(phi), mtot*sinth*sin(phi), mtot*cos(theta) );
   result->SetMomentum(tempVector);

   return result;
}
