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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
//080612 Bug fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #2,3
//080709 Bug fix Sampling Legendre expansion by T. Koi   
//101110 Bug fix in MF=6, LAW=2 case; contribution from E. Mendoza, D. Cano-Ott (CIEMAT)
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPDiscreteTwoBody.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4ParticleHPVector.hh"
#include "G4ParticleHPLegendreStore.hh"

G4ReactionProduct * G4ParticleHPDiscreteTwoBody::Sample(G4double anEnergy, G4double massCode, G4double )
{ // Interpolation still only for the most used parts; rest to be Done @@@@@
   G4ReactionProduct * result = new G4ReactionProduct;
   G4int Z = static_cast<G4int>(massCode/1000);
   G4int A = static_cast<G4int>(massCode-1000*Z);

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
     if(Z!=2) throw G4HadronicException(__FILE__, __LINE__, "Unknown ion case 1");    
   }
   else
   {
     throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPDiscreteTwoBody: Unknown ion case 2");
   }
   
// get cosine(theta)
   G4int i(0), it(0);
   G4double cosTh(0);
   for(i=0; i<nEnergy; i++)
   {
     it = i;
     if(theCoeff[i].GetEnergy()>anEnergy) break;
   }
   if(it==0||it==nEnergy-1)
   {
     if(theCoeff[it].GetRepresentation()==0)
     {
//TK Legendre expansion
       G4ParticleHPLegendreStore theStore(1);
       theStore.SetCoeff(0, theCoeff);
       theStore.SetManager(theManager);
       //cosTh = theStore.SampleMax(anEnergy);
       //080612TK contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #3
       cosTh = theStore.SampleDiscreteTwoBody(anEnergy);
     }
     else if(theCoeff[it].GetRepresentation()==12) // means LINLIN
     {
       G4ParticleHPVector theStore; 
       G4InterpolationManager aManager;
       aManager.Init(LINLIN, theCoeff[it].GetNumberOfPoly()/2);
       theStore.SetInterpolationManager(aManager);
       for(i=0;i<theCoeff[it].GetNumberOfPoly(); i+=2)
       {
         //101110
         //theStore.SetX(i, theCoeff[it].GetCoeff(i));
         //theStore.SetY(i, theCoeff[it].GetCoeff(i));
         theStore.SetX(i/2, theCoeff[it].GetCoeff(i));
         theStore.SetY(i/2, theCoeff[it].GetCoeff(i+1));
       }
       cosTh = theStore.Sample();
     }
     else if(theCoeff[it].GetRepresentation()==14) //this is LOGLIN
     {
       G4ParticleHPVector theStore;
       G4InterpolationManager aManager;
       aManager.Init(LOGLIN, theCoeff[it].GetNumberOfPoly()/2);
       theStore.SetInterpolationManager(aManager);
       for(i=0;i<theCoeff[it].GetNumberOfPoly(); i+=2)
       {
         //101110
         //theStore.SetX(i, theCoeff[it].GetCoeff(i));
         //theStore.SetY(i, theCoeff[it].GetCoeff(i));
         theStore.SetX(i/2, theCoeff[it].GetCoeff(i));
         theStore.SetY(i/2, theCoeff[it].GetCoeff(i+1));
       }
       cosTh = theStore.Sample(); 
     }
     else
     {
       throw G4HadronicException(__FILE__, __LINE__, "unknown representation type in Two-body scattering");
     }
   }
   else
   {
     if(!bCheckDiffCoeffRepr || theCoeff[it].GetRepresentation() == theCoeff[it-1].GetRepresentation())
     {
       if(theCoeff[it].GetRepresentation()==0)
       {
//TK Legendre expansion
	 G4ParticleHPLegendreStore theStore(2);
	 theStore.SetCoeff(0, &(theCoeff[it-1]));
	 theStore.SetCoeff(1, &(theCoeff[it]));
         G4InterpolationManager aManager;
         aManager.Init(theManager.GetScheme(it), 2);
         theStore.SetManager(aManager);
	 //cosTh = theStore.SampleMax(anEnergy);
//080709 TKDB
         cosTh = theStore.SampleDiscreteTwoBody(anEnergy);
       }
       else if(theCoeff[it].GetRepresentation()==12) // LINLIN
       {
	 G4ParticleHPVector theBuff1;
         G4InterpolationManager aManager1;
         aManager1.Init(LINLIN, theCoeff[it-1].GetNumberOfPoly()/2);
         theBuff1.SetInterpolationManager(aManager1);
	 for(i=0;i<theCoeff[it-1].GetNumberOfPoly(); i+=2)
	 {
           //101110
           //theBuff1.SetX(i, theCoeff[it-1].GetCoeff(i));
           //theBuff1.SetY(i, theCoeff[it-1].GetCoeff(i));
           theBuff1.SetX(i/2, theCoeff[it-1].GetCoeff(i));
           theBuff1.SetY(i/2, theCoeff[it-1].GetCoeff(i+1));
	 }
	 G4ParticleHPVector theBuff2;
         G4InterpolationManager aManager2;
         aManager2.Init(LINLIN, theCoeff[it].GetNumberOfPoly()/2);
         theBuff2.SetInterpolationManager(aManager2);
	 for(i=0;i<theCoeff[it].GetNumberOfPoly(); i+=2)
	 {
           //theBuff2.SetX(i, theCoeff[it].GetCoeff(i));
           //theBuff2.SetY(i, theCoeff[it].GetCoeff(i));
           theBuff2.SetX(i/2, theCoeff[it].GetCoeff(i));
           theBuff2.SetY(i/2, theCoeff[it].GetCoeff(i+1));
	 }

	 G4double x1 = theCoeff[it-1].GetEnergy();
	 G4double x2 = theCoeff[it].GetEnergy();
	 G4double x = anEnergy;
	 G4double y1, y2, y, mu;

	 G4ParticleHPVector theStore1;
         theStore1.SetInterpolationManager(aManager1);
	 G4ParticleHPVector theStore2;
         theStore2.SetInterpolationManager(aManager2);
	 G4ParticleHPVector theStore;
	 
	 // for fixed mu get p1, p2 and interpolate according to x
         for(i=0; i<theBuff1.GetVectorLength(); i++)
         {
           mu = theBuff1.GetX(i);
           y1 = theBuff1.GetY(i);
           y2 = theBuff2.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore1.SetData(i, mu, y);
         }
         for(i=0; i<theBuff2.GetVectorLength(); i++)
         {
           mu = theBuff2.GetX(i);
           y1 = theBuff2.GetY(i);
           y2 = theBuff1.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore2.SetData(i, mu, y);
         }
         theStore.Merge(&theStore1, &theStore2); // merge takes care of interpolationschemes
	 cosTh = theStore.Sample();
       }
       else if(theCoeff[it].GetRepresentation()==14) //TK LOG_LIN
       {
	 G4ParticleHPVector theBuff1;
         G4InterpolationManager aManager1;
         aManager1.Init(LOGLIN, theCoeff[it-1].GetNumberOfPoly()/2);
         theBuff1.SetInterpolationManager(aManager1);
	 for(i=0;i<theCoeff[it-1].GetNumberOfPoly(); i+=2)
	 {
           //101110
           //theBuff1.SetX(i, theCoeff[it-1].GetCoeff(i));
           //theBuff1.SetY(i, theCoeff[it-1].GetCoeff(i));
           theBuff1.SetX(i/2, theCoeff[it-1].GetCoeff(i));
           theBuff1.SetY(i/2, theCoeff[it-1].GetCoeff(i+1));
	 }
	 
	 G4ParticleHPVector theBuff2;
         G4InterpolationManager aManager2;
         aManager2.Init(LOGLIN, theCoeff[it].GetNumberOfPoly()/2);
         theBuff2.SetInterpolationManager(aManager2);
	 for(i=0;i<theCoeff[it].GetNumberOfPoly(); i+=2)
	 {
           //101110
           //theBuff2.SetX(i, theCoeff[it].GetCoeff(i));
           //theBuff2.SetY(i, theCoeff[it].GetCoeff(i));
           theBuff2.SetX(i/2, theCoeff[it].GetCoeff(i));
           theBuff2.SetY(i/2, theCoeff[it].GetCoeff(i+1));
	 }

	 G4double x1 = theCoeff[it-1].GetEnergy();
	 G4double x2 = theCoeff[it].GetEnergy();
	 G4double x = anEnergy;
	 G4double y1, y2, y, mu;

	 G4ParticleHPVector theStore1;
         theStore1.SetInterpolationManager(aManager1);
	 G4ParticleHPVector theStore2;
         theStore2.SetInterpolationManager(aManager2);
	 G4ParticleHPVector theStore;
	 
	 // for fixed mu get p1, p2 and interpolate according to x
         for(i=0; i<theBuff1.GetVectorLength(); i++)
         {
           mu = theBuff1.GetX(i);
           y1 = theBuff1.GetY(i);
           y2 = theBuff2.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore1.SetData(i, mu, y);
         }
         for(i=0; i<theBuff2.GetVectorLength(); i++)
         {
           mu = theBuff2.GetX(i);
           y1 = theBuff2.GetY(i);
           y2 = theBuff1.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore2.SetData(i, mu, y);
         }
         theStore.Merge(&theStore1, &theStore2); 
	 cosTh = theStore.Sample(); 
       }
       else
       {
         throw G4HadronicException(__FILE__, __LINE__, "Two neighbouring distributions with different interpolation");
       }
     }
     else
     {
       G4cout << " theCoeff[it].GetRepresent MEM " << it << " " << &theCoeff[it] <<  "  " << &theCoeff[it-1] << G4endl;
       G4cout << " theCoeff[it].GetRepresent " << it << " " << theCoeff[it].GetRepresentation() <<  " != " << theCoeff[it-1].GetRepresentation() << G4endl;

       throw G4HadronicException(__FILE__, __LINE__, "unknown representation type in Two-body scattering, case 2");
     }
   }
   
// now get the energy from kinematics and Q-value.

   //G4double restEnergy = anEnergy+GetQValue();
   
// assumed to be in CMS @@@@@@@@@@@@@@@@@

   //080612TK contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #2
   //G4double residualMass =   GetTarget()->GetMass() + GetNeutron()->GetMass()
   //                        - result->GetMass() - GetQValue();
   //G4double kinE = restEnergy/(1+result->GetMass()/residualMass); // non relativistic @@
   G4double A1     =  GetTarget()->GetMass()/GetProjectileRP()->GetMass(); 
   G4double A1prim =  result->GetMass()/GetProjectileRP()->GetMass();
   //G4double E1     =  (A1+1)*(A1+1)/A1/A1*anEnergy; 
   //Bug fix Bugzilla #1815
   G4double E1     =  anEnergy; 
   G4double kinE = (A1+1-A1prim)/(A1+1)/(A1+1)*(A1*E1+(1+A1)*GetQValue());

   result->SetKineticEnergy(kinE); // non relativistic @@
   G4double phi = CLHEP::twopi*G4UniformRand();
   G4double theta = std::acos(cosTh);
   G4double sinth = std::sin(theta);
   G4double mtot = result->GetTotalMomentum(); 
   G4ThreeVector tempVector(mtot*sinth*std::cos(phi), mtot*sinth*std::sin(phi), mtot*std::cos(theta) );
   result->SetMomentum(tempVector);
   
// some garbage collection
   
// return the result   
   return result;
}
