// This code implementation is the intellectual property of
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPNBodyPhaseSpace.cc,v 1.2 1999-07-02 10:01:14 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4NeutronHPNBodyPhaseSpace.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

G4ReactionProduct * G4NeutronHPNBodyPhaseSpace::Sample(G4double anEnergy, G4double massCode, G4double mass)
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
     G4Exception("G4NeutronHPNBodyPhaseSpace: Unknown ion case 2");
   }

// Get the energy from phase-space distribution
   // in CMS
   // P = Cn*sqrt(E')*(Emax-E')**(3*n/2-4)
   G4double maxE = GetEmax(anEnergy, result->GetMass());
   G4double energy;
   G4double max;
   G4double tMass = GetTarget()->GetMass();
   G4double pMass = GetNeutron()->GetMass();
   if(theTotalCount<=3)
   {
     max = maxE/2.;
   }
   else if(theTotalCount==4)
   {
     max = maxE/5.;
   }
   else if(theTotalCount==5)
   {
     max = maxE/8.;
   }
   else
   {
     G4Exception("NeutronHP Phase-space distribution cannot cope with this number of particles");
   }
   G4double testit;
   G4double rand0 = Prob(max, maxE, theTotalCount);
   G4double rand;
   
   do
   {
     rand = rand0*G4UniformRand();
     energy = maxE*G4UniformRand();
     testit = Prob(energy, maxE, theTotalCount);
   }
   while(rand > testit);
   result->SetKineticEnergy(energy);
   
// now do random direction
   G4double cosTh = 2.*G4UniformRand()-1.;
   G4double phi = twopi*G4UniformRand();
   G4double theta = acos(cosTh);
   G4double sinth = sin(theta);
   G4double mtot = result->GetTotalMomentum(); 
   G4ThreeVector tempVector(mtot*sinth*cos(phi), mtot*sinth*sin(phi), mtot*cos(theta) );
   result->SetMomentum(tempVector);
   G4ReactionProduct aCMS = *GetTarget()+*GetNeutron();
   result->Lorentz(*result, -1.*aCMS);
   return result;
}
